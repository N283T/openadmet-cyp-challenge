"""Summarize held-out D6 errors without selecting a model or validation gate."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def summarize(frame: pl.DataFrame, population: str, fold: int, report_st_rae: bool) -> dict:
    y = frame["y_true"].to_numpy()
    pred = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    soft = np.maximum(pred - hi, 0) + np.maximum(lo - pred, 0)
    baseline = np.maximum(y.mean() - hi, 0) + np.maximum(lo - y.mean(), 0)
    denominator = float(baseline.sum())
    return {
        "population": population,
        "fold": fold,
        "n": len(y),
        "mean_pic50": float(y.mean()),
        "median_pic50": float(np.median(y)),
        "mean_prediction": float(pred.mean()),
        "mean_signed_error": float((pred - y).mean()),
        "mean_soft_error": float(soft.mean()),
        "soft_error_sum": float(soft.sum()),
        "soft_baseline_sum": denominator,
        "ST_RAE": (
            float(soft.sum() / denominator)
            if report_st_rae and denominator > 0 else None
        ),
        "MAE": float(np.abs(y - pred).mean()),
        "Spearman": float(spearmanr(y, pred).statistic) if len(y) >= 3 else None,
    }


def main() -> None:
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = config["endpoint"]
    keys = ["compound_id", "endpoint", "fold"]
    reference = pl.read_parquet(inputs / "reference-oof.parquet").filter(
        pl.col("endpoint") == endpoint
    ).select(*keys, "y_true", "conf_low", "conf_high", "is_neighbor", "split_id")
    control = pl.read_parquet(inputs / "control-oof.parquet").filter(
        (pl.col("endpoint") == endpoint) & (pl.col("scenario") == config["control_scenario"])
    ).select(*keys, "prediction")
    if reference.height != 1457 or control.height != reference.height:
        raise ValueError("Unexpected current D6 OOF coverage")
    if reference.select(keys).is_duplicated().any() or control.select(keys).is_duplicated().any():
        raise ValueError("Duplicate OOF keys")
    pairs = pl.read_parquet(inputs / "analog-pairs.parquet")
    analogs = pairs.select(pl.col("analog_id").alias("compound_id")).unique()
    aligned = reference.join(control, on=keys, validate="1:1").join(
        analogs.with_columns(pl.lit(True).alias("in_proxy")),
        on="compound_id", how="left", validate="m:1",
    ).with_columns(pl.col("in_proxy").fill_null(False))
    if aligned.height != 1457 or aligned["split_id"].n_unique() != 1:
        raise ValueError("OOF alignment or split mismatch")
    low = config["low_activity_cutoff"]
    high = config["high_activity_cutoff"]
    subsets = {
        "all": aligned,
        "other_cyp_analog_proxy": aligned.filter(pl.col("in_proxy")),
        "outside_proxy": aligned.filter(~pl.col("in_proxy")),
        "below_4": aligned.filter(pl.col("y_true") < low),
        "4_to_5": aligned.filter((pl.col("y_true") >= low) & (pl.col("y_true") < high)),
        "at_least_5": aligned.filter(pl.col("y_true") >= high),
        "proxy_below_4": aligned.filter(pl.col("in_proxy") & (pl.col("y_true") < low)),
    }
    rows = []
    for name, frame in subsets.items():
        if frame.is_empty():
            continue
        report_st_rae = name in {"all", "other_cyp_analog_proxy", "outside_proxy"}
        rows.append(summarize(frame, name, -1, report_st_rae))
        for fold in range(5):
            part = frame.filter(pl.col("fold") == fold)
            if not part.is_empty():
                rows.append(summarize(part, name, fold, report_st_rae))
    pl.DataFrame(rows).write_csv(outputs / "metrics.csv")
    aligned.write_parquet(outputs / "aligned-oof.parquet")
    (outputs / "summary.json").write_text(json.dumps({
        "purpose": "Descriptive only; proxy and outcome strata cannot select a candidate",
        "st_rae_policy": (
            "Omit ST-RAE on outcome-conditioned activity strata, where the "
            "mean baseline denominator can collapse"
        ),
        "oof_rows": aligned.height,
        "proxy_rows": subsets["other_cyp_analog_proxy"].height,
        "split_id": aligned["split_id"][0],
        "config": config,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
