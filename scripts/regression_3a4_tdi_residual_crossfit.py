"""Test whether held-out classification A4 scores add residual regression signal."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


def st_rae(y: np.ndarray, pred: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> float:
    center = y.mean()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float((np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum() / denominator)


def metrics(frame: pl.DataFrame, pred: np.ndarray) -> dict[str, float]:
    y = frame["y_true"].to_numpy()
    return {
        "ST_RAE": st_rae(y, pred, frame["conf_low"].to_numpy(), frame["conf_high"].to_numpy()),
        "MAE": float(np.abs(y - pred).mean()),
        "Spearman": float(spearmanr(y, pred).statistic),
    }


def parity(group: str, repeat: int) -> int:
    digest = hashlib.sha256(f"{repeat}:{group}".encode()).digest()
    return digest[0] & 1


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    reg = (
        pl.read_parquet(inp / "regression-oof.parquet")
        .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["base_method"]))
        .sort("compound_id")
    )
    cls = (
        pl.read_parquet(inp / "classification-oof.parquet")
        .filter(pl.col("model") == cfg["tdi_model"])
        .select("compound_id", pl.col("prediction").alias("tdi_probability"))
    )
    assignments = pl.read_parquet(inp / "parent-assignments.parquet").select(
        "compound_id", "group", pl.col("fold").alias("assigned_fold"), "is_anchor"
    )
    data = (
        reg.join(cls, on="compound_id", validate="1:1")
        .join(assignments, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    if data.height != 2300 or data["compound_id"].n_unique() != 2300:
        raise ValueError("Incomplete 3A4 regression/classification OOF overlap")
    if not (data["fold"] == data["assigned_fold"]).all() or data["is_anchor"].any():
        raise ValueError("Regression held-fold membership mismatch")
    if not np.isfinite(data.select("prediction", "tdi_probability").to_numpy()).all():
        raise ValueError("Nonfinite OOF predictors")

    y = data["y_true"].to_numpy()
    base = data["prediction"].to_numpy()
    tdi = data["tdi_probability"].to_numpy()
    fold = data["fold"].to_numpy()
    groups = data["group"].to_list()
    methods = ["base_only", "base_plus_tdi"]
    predictions = []
    score_rows = []
    fit_rows = []
    for repeat in range(cfg["repeats"]):
        half = np.asarray([parity(group, repeat) for group in groups])
        corrected = {method: np.full(data.height, np.nan) for method in methods}
        for held_fold in range(5):
            for held_half in [0, 1]:
                train = (fold == held_fold) & (half != held_half)
                valid = (fold == held_fold) & (half == held_half)
                if train.sum() < 100 or valid.sum() < 100:
                    raise ValueError("Insufficient disjoint within-fold calibration halves")
                for method in methods:
                    x = base[:, None] if method == "base_only" else np.column_stack([base, tdi])
                    scaler = StandardScaler().fit(x[train])
                    model = Ridge(alpha=cfg["ridge_alpha"]).fit(
                        scaler.transform(x[train]), (y - base)[train]
                    )
                    corrected[method][valid] = model.predict(scaler.transform(x[valid]))
                    fit_rows.append(
                        {
                            "repeat": repeat,
                            "fold": held_fold,
                            "calibration_half": 1 - held_half,
                            "method": method,
                            "train": int(train.sum()),
                            "valid": int(valid.sum()),
                            "intercept": float(model.intercept_),
                            "coefficients": json.dumps(model.coef_.tolist()),
                        }
                    )
        for method, residual in corrected.items():
            if not np.isfinite(residual).all():
                raise ValueError(f"Incomplete residual prediction: {method}")
            for strength in cfg["correction_strengths"]:
                name = f"{method}_{strength:g}"
                prediction = base + strength * residual
                predictions.append(
                    data.select(
                        "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                    ).with_columns(
                        pl.lit(repeat).alias("repeat"),
                        pl.lit(name).alias("method"),
                        pl.Series("prediction", prediction),
                    )
                )
                score_rows.append(
                    {
                        "repeat": repeat,
                        "method": name,
                        "n": data.height,
                        **metrics(data, prediction),
                    }
                )
        score_rows.append(
            {"repeat": repeat, "method": "identity", "n": data.height, **metrics(data, base)}
        )
    pl.concat(predictions).write_parquet(out / "crossfit-oof.parquet")
    pl.DataFrame(score_rows).sort("repeat", "method").write_csv(out / "repeat-metrics.csv")
    pl.DataFrame(fit_rows).sort("repeat", "fold", "calibration_half", "method").write_csv(
        out / "fits.csv"
    )
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "rows": data.height,
                "regression_folds": 5,
                "group_disjoint_within_fold": True,
                "tdi_probability_is_classification_oof": True,
                "classification_regression_fold_maps_differ": True,
                "qualification": "Within-fold calibration/evaluation labels are disjoint. "
                "OOF classifier scores exclude each compound's own TDI label. "
                "Multiple repeats and fixed models remain exploratory; test transfer unverified.",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
