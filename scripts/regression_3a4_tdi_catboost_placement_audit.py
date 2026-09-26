"""Check whether the small 3A4 placed blend gain persists by fold and series."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def st_rae(lo: np.ndarray, hi: np.ndarray, y: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold"]
    base = (
        pl.read_parquet(inp / "scored-bank-oof.parquet")
        .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["base_method"]))
        .select(keys + ["y_true", "conf_low", "conf_high", "prediction"])
    )
    placed = pl.read_parquet(inp / "placed-oof.parquet").select(keys + ["placed_prediction"])
    groups = pl.read_parquet(inp / "parent-assignments.parquet").select(["compound_id", "group"])
    data = (
        base.join(placed, on=keys, how="inner", validate="1:1")
        .join(groups, on="compound_id", how="inner", validate="1:1")
        .sort(keys)
    )
    if data.height != 2300 or data.null_count().row(0) != tuple(0 for _ in data.columns):
        raise ValueError("Incomplete 3A4 OOF or group alignment")
    y, lo, hi = (data[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"])
    old, new = (data[c].to_numpy() for c in ["prediction", "placed_prediction"])
    folds = data["fold"].to_numpy()
    group_ids = data["group"].to_numpy()
    unique_groups, inverse = np.unique(group_ids, return_inverse=True)
    group_rows = [np.flatnonzero(inverse == i) for i in range(len(unique_groups))]
    rng = np.random.default_rng(cfg["seed"])
    fold_rows = []
    bootstrap_rows = []
    for weight in cfg["weights"]:
        blended = (1 - weight) * old + weight * new
        for fold in np.unique(folds):
            mask = folds == fold
            fold_rows.append(
                {
                    "weight": weight,
                    "fold": int(fold),
                    "n": int(mask.sum()),
                    "base_ST_RAE": st_rae(lo[mask], hi[mask], y[mask], old[mask]),
                    "blend_ST_RAE": st_rae(lo[mask], hi[mask], y[mask], blended[mask]),
                }
            )
        for repeat in range(cfg["bootstrap_repeats"]):
            sampled = rng.integers(0, len(unique_groups), len(unique_groups))
            idx = np.concatenate([group_rows[i] for i in sampled])
            delta = st_rae(lo[idx], hi[idx], y[idx], blended[idx]) - st_rae(
                lo[idx], hi[idx], y[idx], old[idx]
            )
            bootstrap_rows.append({"weight": weight, "repeat": repeat, "delta_ST_RAE": delta})
    folds_df = pl.DataFrame(fold_rows).with_columns(
        (pl.col("blend_ST_RAE") - pl.col("base_ST_RAE")).alias("delta_ST_RAE")
    )
    boot_df = pl.DataFrame(bootstrap_rows)
    folds_df.write_csv(out / "fold-metrics.csv")
    boot_df.write_csv(out / "bootstrap.csv")
    summary = {
        "n": data.height,
        "n_groups": len(unique_groups),
        "bootstrap_repeats": cfg["bootstrap_repeats"],
        "weights": {},
    }
    for weight in cfg["weights"]:
        values = boot_df.filter(pl.col("weight") == weight)["delta_ST_RAE"].to_numpy()
        overall = (1 - weight) * old + weight * new
        summary["weights"][str(weight)] = {
            "overall_delta_ST_RAE": st_rae(lo, hi, y, overall) - st_rae(lo, hi, y, old),
            "folds_improved": int(
                folds_df.filter((pl.col("weight") == weight) & (pl.col("delta_ST_RAE") < 0)).height
            ),
            "bootstrap_delta_95pct": np.quantile(values, [0.025, 0.975]).tolist(),
            "bootstrap_fraction_improved": float(np.mean(values < 0)),
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
