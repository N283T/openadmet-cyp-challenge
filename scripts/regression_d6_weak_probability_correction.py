"""Test fixed, outer-training-centered weak activity probability corrections."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    error = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    mean = y.mean()
    baseline = (np.maximum(lo - mean, 0) + np.maximum(mean - hi, 0)).sum()
    return float(error / baseline)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    classifier = pl.read_parquet(inp / "classifier-oof.parquet").filter(
        (pl.col("endpoint") == endpoint) & (pl.col("view") == cfg["classification_view"])
        & pl.col("label_rule").is_in(cfg["label_rules"])
    )
    labels = pl.read_parquet(inp / "measurements.parquet").filter(
        pl.col("endpoint") == endpoint
    ).select("compound_id", pl.col("value").alias("y_true"), "conf_high")
    if labels.height != 1493:
        raise ValueError("Unexpected D6 labels")
    results = []
    fold_params = []
    for split in ["parent", "scaffold"]:
        base = pl.read_parquet(inp / f"{split}-base-oof.parquet").filter(
            pl.col("endpoint") == endpoint
        ).select("compound_id", "endpoint", "fold", "y_true", "conf_low",
                 "conf_high", "prediction")
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet")
        training = labels.join(assignments.select(
            "compound_id", "fold", "is_anchor"
        ), on="compound_id", validate="1:1")
        if base.height != 1457 or training.height != 1493:
            raise ValueError(f"Unexpected {split} coverage")
        if base.select("compound_id", "endpoint", "fold").is_duplicated().any():
            raise ValueError(f"Duplicate {split} OOF keys")
        for rule in cfg["label_rules"]:
            probabilities = classifier.filter(
                (pl.col("split_id") == split) & (pl.col("label_rule") == rule)
            ).select("compound_id", "endpoint", "fold", "probability")
            joined = base.join(probabilities, on=["compound_id", "endpoint", "fold"],
                               validate="1:1").sort("compound_id")
            if joined.height != 1457:
                raise ValueError(f"Incomplete {split}/{rule} probability alignment")
            if not np.allclose(joined["y_true"], classifier.filter(
                (pl.col("split_id") == split) & (pl.col("label_rule") == rule)
            ).sort("compound_id")["y_true"]):
                raise ValueError("Classifier/regressor label mismatch")
            for fold in range(5):
                held = joined.filter(pl.col("fold") == fold)
                fit = training.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
                weak = (
                    fit["y_true"].to_numpy() if rule == "point"
                    else fit["conf_high"].to_numpy()
                ) < cfg["weak_cutoff"]
                y_fit = fit["y_true"].to_numpy()
                prevalence = float(weak.mean())
                delta = float(y_fit[weak].mean() - y_fit[~weak].mean())
                fold_params.append({
                    "split_id": split, "label_rule": rule, "fold": fold,
                    "train_n": fit.height, "weak_n": int(weak.sum()),
                    "prevalence": prevalence, "mean_gap": delta,
                })
                signal = (held["probability"].to_numpy() - prevalence) * delta
                for name, strength in [("base", 0.0)] + [
                    (f"strength_{s:g}", s) for s in cfg["strengths"]
                ]:
                    results.append(held.select(
                        "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                    ).with_columns(
                        pl.lit(split).alias("split_id"),
                        pl.lit(rule).alias("label_rule"),
                        pl.lit(name).alias("method"),
                        pl.Series("prediction", held["prediction"].to_numpy() + strength * signal),
                    ))
    oof = pl.concat(results).sort("split_id", "label_rule", "method", "compound_id")
    expected = 2 * len(cfg["label_rules"]) * (1 + len(cfg["strengths"])) * 1457
    if oof.height != expected or oof.select(
        "split_id", "label_rule", "method", "compound_id"
    ).is_duplicated().any():
        raise ValueError("Incomplete correction OOF")
    metrics = []
    for (split, rule, method), frame in oof.group_by("split_id", "label_rule", "method"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            part = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y = part["y_true"].to_numpy()
            lo = part["conf_low"].to_numpy()
            hi = part["conf_high"].to_numpy()
            pred = part["prediction"].to_numpy()
            low = y < cfg["weak_cutoff"]
            metrics.append({
                "split_id": split, "label_rule": rule, "method": method,
                "fold": fold, "n": part.height,
                "ST_RAE": st_rae(y, lo, hi, pred),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
                "weak_n": int(low.sum()),
                "weak_MAE": float(np.abs(y[low] - pred[low]).mean()),
                "weak_mean_soft_error": float((
                    np.maximum(lo[low] - pred[low], 0)
                    + np.maximum(pred[low] - hi[low], 0)
                ).mean()),
            })
    pl.DataFrame(metrics).sort("split_id", "label_rule", "method", "fold").write_csv(
        out / "metrics.csv"
    )
    pl.DataFrame(fold_params).write_csv(out / "fold-parameters.csv")
    oof.write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg, "rows": oof.height,
        "qualification": "Exploratory post-classifier screen; no independent confirmation",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
