"""Fit identical CatBoost readouts with prespecified weak-label weights."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostRegressor
from scipy.stats import spearmanr


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    center = y.mean()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / "measurements.parquet").filter(
        pl.col("endpoint") == cfg["endpoint"]
    ).select("compound_id", "endpoint", pl.col("value").alias("y_true"), "conf_low", "conf_high")
    if labels.height != 1493 or labels["compound_id"].n_unique() != labels.height:
        raise ValueError("Unexpected D6 label coverage")
    features = pl.read_parquet(inp / "scores.parquet")
    columns = sorted(c for c in features.columns if c.startswith("score_"))
    if len(columns) != 47 or features["compound_id"].n_unique() != features.height:
        raise ValueError("Unexpected frozen score features")
    outputs = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet")
        if assignments["compound_id"].n_unique() != assignments.height:
            raise ValueError(f"Duplicate {split} assignment")
        data = labels.join(
            assignments.select("compound_id", "fold", "group", "is_anchor"),
            on="compound_id", validate="1:1",
        ).join(features.select("compound_id", *columns), on="compound_id", validate="1:1")
        if data.height != 1493 or data.filter(~pl.col("is_anchor")).height != 1457:
            raise ValueError(f"Incomplete {split} feature alignment")
        data = data.sort("compound_id")
        matrix = data.select(columns).to_numpy().astype(np.float32)
        if not np.isfinite(matrix).all():
            raise ValueError(f"Nonfinite {split} features")
        y = data["y_true"].to_numpy()
        lo = data["conf_low"].to_numpy()
        hi = data["conf_high"].to_numpy()
        folds = data["fold"].to_numpy()
        anchor = data["is_anchor"].to_numpy()
        for fold in range(5):
            train = (folds != fold) | anchor
            valid = (folds == fold) & ~anchor
            if not valid.any() or np.any(train & valid):
                raise ValueError(f"Invalid {split} fold {fold}")
            for recipe in cfg["recipes"]:
                weights = np.ones(int(train.sum()), dtype=float)
                if recipe["rule"] == "point":
                    weights[y[train] < cfg["weak_cutoff"]] = recipe["weight"]
                elif recipe["rule"] == "upper":
                    weights[hi[train] < cfg["weak_cutoff"]] = recipe["weight"]
                elif recipe["rule"] != "none":
                    raise ValueError(recipe["rule"])
                settings = cfg["model"]
                model = CatBoostRegressor(
                    iterations=settings["iterations"], depth=settings["depth"],
                    learning_rate=settings["learning_rate"],
                    l2_leaf_reg=settings["l2_leaf_reg"],
                    loss_function=settings["loss_function"],
                    random_seed=settings["seed"] + fold,
                    thread_count=settings["thread_count"], verbose=False,
                    allow_writing_files=False,
                )
                model.fit(matrix[train], y[train], sample_weight=weights)
                pred = model.predict(matrix[valid]).astype(float)
                rows = data.filter(pl.Series(valid)).select(
                    "compound_id", "endpoint", "fold", "group", "y_true",
                    "conf_low", "conf_high",
                ).with_columns(
                    pl.lit(split).alias("split_id"),
                    pl.lit(recipe["name"]).alias("method"),
                    pl.Series("prediction", pred),
                )
                outputs.append(rows)
                print(split, fold, recipe["name"], "n", len(pred), flush=True)
    oof = pl.concat(outputs).sort("split_id", "method", "compound_id")
    expected = len(cfg["splits"]) * len(cfg["recipes"]) * 1457
    duplicate = oof.select("split_id", "method", "compound_id").is_duplicated().any()
    if oof.height != expected or duplicate:
        raise ValueError("Incomplete or duplicate OOF output")
    metric_rows = []
    for (split, method), frame in oof.group_by("split_id", "method"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            part = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y = part["y_true"].to_numpy()
            lo = part["conf_low"].to_numpy()
            hi = part["conf_high"].to_numpy()
            pred = part["prediction"].to_numpy()
            low = y < cfg["weak_cutoff"]
            metric_rows.append({
                "split_id": split, "method": method, "fold": fold, "n": len(y),
                "ST_RAE": st_rae(y, lo, hi, pred),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
                "weak_n": int(low.sum()),
                "weak_MAE": float(np.abs(y[low] - pred[low]).mean()),
                "weak_mean_signed_error": float((pred[low] - y[low]).mean()),
                "weak_mean_soft_error": float((
                    np.maximum(lo[low] - pred[low], 0)
                    + np.maximum(pred[low] - hi[low], 0)
                ).mean()),
            })
    pl.DataFrame(metric_rows).sort("split_id", "method", "fold").write_csv(out / "metrics.csv")
    oof.write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg,
        "rule": "Training weights use outer-training labels only; activity strata are diagnostics",
        "rows": oof.height,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
