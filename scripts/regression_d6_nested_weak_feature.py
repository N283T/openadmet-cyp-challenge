"""Compare matched regressors with and without an inner-cross-fitted weak signal."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostClassifier, CatBoostRegressor
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    error = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    center = y.mean()
    baseline = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(error / baseline)


def classifier(config: dict, seed: int, threads: int) -> CatBoostClassifier:
    return CatBoostClassifier(
        iterations=config["iterations"], depth=config["depth"],
        learning_rate=config["learning_rate"], l2_leaf_reg=config["l2_leaf_reg"],
        loss_function="Logloss", random_seed=seed, thread_count=threads,
        verbose=False, allow_writing_files=False,
    )


def regressor(config: dict, seed: int, threads: int) -> CatBoostRegressor:
    return CatBoostRegressor(
        iterations=config["iterations"], depth=config["depth"],
        learning_rate=config["learning_rate"], l2_leaf_reg=config["l2_leaf_reg"],
        loss_function=config["loss_function"], random_seed=seed,
        thread_count=threads, verbose=False, allow_writing_files=False,
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / "measurements.parquet").filter(
        pl.col("endpoint") == cfg["endpoint"]
    ).select("compound_id", "endpoint", pl.col("value").alias("y_true"),
             "conf_low", "conf_high")
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    columns = sorted(c for c in scores.columns if c.startswith("score_"))
    emb_columns = sorted(c for c in medium.columns if c.startswith("emb_"))
    if labels.height != 1493 or len(columns) != 47 or len(emb_columns) != 600:
        raise ValueError("Unexpected D6 inputs")
    features = scores.select("compound_id", *columns).join(
        medium.select("compound_id", *emb_columns), on="compound_id", validate="1:1"
    )
    columns += emb_columns
    outputs = []
    fit_rows = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet")
        data = labels.join(assignments.select(
            "compound_id", "fold", "group", "is_anchor"
        ), on="compound_id", validate="1:1").join(
            features, on="compound_id", validate="1:1"
        ).sort("compound_id")
        if data.height != 1493 or data.filter(~pl.col("is_anchor")).height != 1457:
            raise ValueError(f"Incomplete {split} coverage")
        x = data.select(columns).to_numpy().astype(np.float32)
        if not np.isfinite(x).all():
            raise ValueError(f"Nonfinite {split} features")
        y = data["y_true"].to_numpy()
        weak = (y < cfg["weak_cutoff"]).astype(int)
        fold_ids = data["fold"].to_numpy()
        anchor = data["is_anchor"].to_numpy()
        groups = data["group"].to_numpy()
        for outer in range(5):
            train_idx = np.flatnonzero((fold_ids != outer) | anchor)
            valid_idx = np.flatnonzero((fold_ids == outer) & ~anchor)
            if np.intersect1d(train_idx, valid_idx).size:
                raise ValueError("Outer train/validation overlap")
            inner_probability = np.full(len(train_idx), np.nan)
            inner = GroupKFold(n_splits=cfg["inner_group_folds"])
            for inner_fold, (fit_idx, held_idx) in enumerate(inner.split(
                x[train_idx], weak[train_idx], groups[train_idx]
            )):
                if set(groups[train_idx[fit_idx]]) & set(groups[train_idx[held_idx]]):
                    raise ValueError("Inner group leakage")
                model = classifier(cfg["classifier"], cfg["seed"] + outer * 10 + inner_fold,
                                   cfg["threads"])
                model.fit(x[train_idx[fit_idx]], weak[train_idx[fit_idx]])
                inner_probability[held_idx] = model.predict_proba(x[train_idx[held_idx]])[:, 1]
            if not np.isfinite(inner_probability).all():
                raise ValueError("Incomplete inner probabilities")
            full_classifier = classifier(cfg["classifier"], cfg["seed"] + outer,
                                         cfg["threads"])
            full_classifier.fit(x[train_idx], weak[train_idx])
            outer_probability = full_classifier.predict_proba(x[valid_idx])[:, 1]
            for method in ["control", "nested_weak_feature"]:
                if method == "control":
                    train_x, valid_x = x[train_idx], x[valid_idx]
                else:
                    train_x = np.column_stack([x[train_idx], inner_probability])
                    valid_x = np.column_stack([x[valid_idx], outer_probability])
                model = regressor(cfg["regressor"], cfg["seed"] + outer, cfg["threads"])
                model.fit(train_x, y[train_idx])
                prediction = model.predict(valid_x).astype(float)
                outputs.append(data[valid_idx].select(
                    "compound_id", "endpoint", "fold", "group", "y_true",
                    "conf_low", "conf_high"
                ).with_columns(
                    pl.lit(split).alias("split_id"),
                    pl.lit(method).alias("method"),
                    pl.Series("prediction", prediction),
                ))
            fit_rows.append({
                "split_id": split, "fold": outer, "train_n": len(train_idx),
                "valid_n": len(valid_idx), "train_weak": int(weak[train_idx].sum()),
                "valid_weak": int(weak[valid_idx].sum()),
                "mean_inner_probability": float(inner_probability.mean()),
                "mean_outer_probability": float(outer_probability.mean()),
            })
            print(split, outer, flush=True)
    oof = pl.concat(outputs).sort("split_id", "method", "compound_id")
    if oof.height != len(cfg["splits"]) * 2 * 1457 or oof.select(
        "split_id", "method", "compound_id"
    ).is_duplicated().any():
        raise ValueError("Incomplete nested-feature OOF")
    metric_rows = []
    for (split, method), frame in oof.group_by("split_id", "method"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            part = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y = part["y_true"].to_numpy()
            lo = part["conf_low"].to_numpy()
            hi = part["conf_high"].to_numpy()
            pred = part["prediction"].to_numpy()
            weak = y < cfg["weak_cutoff"]
            metric_rows.append({
                "split_id": split, "method": method, "fold": fold, "n": len(y),
                "ST_RAE": st_rae(y, lo, hi, pred),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
                "weak_n": int(weak.sum()),
                "weak_MAE": float(np.abs(y[weak] - pred[weak]).mean()),
                "weak_mean_soft_error": float((
                    np.maximum(lo[weak] - pred[weak], 0)
                    + np.maximum(pred[weak] - hi[weak], 0)
                ).mean()),
            })
    pl.DataFrame(metric_rows).sort("split_id", "method", "fold").write_csv(out / "metrics.csv")
    pl.DataFrame(fit_rows).write_csv(out / "fit-summary.csv")
    oof.write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg, "rows": oof.height,
        "rule": "Classifier feature is group-cross-fitted on each outer-training fold",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
