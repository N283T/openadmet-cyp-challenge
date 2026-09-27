"""Inspect affine placement on existing test-selection-like analog OOF rows."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    analogs = (
        pl.read_parquet(inp / "analog-pairs.parquet")
        .select(pl.col("analog_id").alias("compound_id"))
        .unique()
        .with_columns(pl.lit(True).alias("analog_proxy"))
    )
    assignments = pl.read_parquet(inp / "parent-assignments.parquet").select("compound_id", "group")
    metrics = []
    fold_metrics = []
    bootstraps = []
    summary = {"config": cfg, "cohorts": {}}
    rng = np.random.default_rng(cfg["seed"])
    for cohort in cfg["cohorts"]:
        name = cohort["name"]
        file = inp / f"{name.lower()}-oof.parquet"
        frame = (
            pl.read_parquet(file)
            .filter(
                (pl.col("endpoint") == cohort["endpoint"]) & (pl.col("method") == cohort["method"])
            )
            .select(
                "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high", "prediction"
            )
            .join(analogs, on="compound_id", how="left", validate="1:1")
            .join(assignments, on="compound_id", how="inner", validate="1:1")
            .with_columns(pl.col("analog_proxy").fill_null(False))
            .sort("compound_id")
        )
        if frame.height not in {1372, 2300} or frame.null_count().row(0) != tuple(
            0 for _ in frame.columns
        ):
            raise ValueError(f"Incomplete {name} OOF/group alignment")
        raw = frame["prediction"].to_numpy()
        global_center = float(raw.mean())
        proxy = frame["analog_proxy"].to_numpy()
        summary["cohorts"][name] = {
            "all_n": frame.height,
            "analog_n": int(proxy.sum()),
            "all_mean_label": float(frame["y_true"].mean()),
            "analog_mean_label": float(frame.filter(pl.col("analog_proxy"))["y_true"].mean()),
            "qualification": "Analogs were selected from 1A2/2C9/3A4 hit anchors; diagnostic only.",
        }
        if proxy.sum() < 100:
            raise ValueError(f"Too few analog labels for {name}")
        for population, selected in [("all", np.ones(frame.height, dtype=bool)), ("analog", proxy)]:
            part = frame.filter(pl.Series("selected", selected))
            y = part["y_true"].to_numpy()
            lo = part["conf_low"].to_numpy()
            hi = part["conf_high"].to_numpy()
            folds = part["fold"].to_numpy()
            group_ids = part["group"].to_numpy()
            unique_groups, inverse = np.unique(group_ids, return_inverse=True)
            group_rows = [np.flatnonzero(inverse == i) for i in range(len(unique_groups))]
            draws = [
                np.concatenate(
                    [group_rows[i] for i in rng.integers(0, len(unique_groups), len(unique_groups))]
                )
                for _ in range(cfg["bootstrap_repeats"])
            ]
            reference = raw[selected]
            baseline = st_rae(y, lo, hi, reference)
            for variant in cohort["variants"]:
                pred = (
                    global_center
                    + variant["shift"]
                    + variant["scale"] * (reference - global_center)
                )
                score = st_rae(y, lo, hi, pred)
                metrics.append(
                    {
                        "cohort": name,
                        "population": population,
                        "variant": variant["name"],
                        "n": len(y),
                        "groups": len(unique_groups),
                        "mean_label": float(y.mean()),
                        "baseline_ST_RAE": baseline,
                        "ST_RAE": score,
                        "delta_ST_RAE": score - baseline,
                        "MAE": float(np.abs(y - pred).mean()),
                    }
                )
                for fold in np.unique(folds):
                    mask = folds == fold
                    fold_metrics.append(
                        {
                            "cohort": name,
                            "population": population,
                            "variant": variant["name"],
                            "fold": int(fold),
                            "n": int(mask.sum()),
                            "ST_RAE": st_rae(y[mask], lo[mask], hi[mask], pred[mask]),
                        }
                    )
                if variant["name"] == "control":
                    continue
                for repeat, idx in enumerate(draws):
                    delta = st_rae(y[idx], lo[idx], hi[idx], pred[idx]) - st_rae(
                        y[idx], lo[idx], hi[idx], reference[idx]
                    )
                    bootstraps.append(
                        {
                            "cohort": name,
                            "population": population,
                            "variant": variant["name"],
                            "repeat": repeat,
                            "delta_ST_RAE": delta,
                        }
                    )
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    pl.DataFrame(fold_metrics).write_csv(out / "fold-metrics.csv")
    pl.DataFrame(bootstraps).write_csv(out / "bootstrap.csv")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
