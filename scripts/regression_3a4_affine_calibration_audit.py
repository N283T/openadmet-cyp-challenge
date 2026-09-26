"""Audit fixed 3A4 affine corrections by fold and molecular series."""

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
    base = (
        pl.read_parquet(inp / "scored-bank-oof.parquet")
        .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["oof_method"]))
        .select("compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high", "prediction")
    )
    group = pl.read_parquet(inp / "parent-assignments.parquet").select("compound_id", "group")
    data = base.join(group, on="compound_id", how="inner", validate="1:1").sort("compound_id")
    if data.height != 2300 or data.null_count().row(0) != tuple(0 for _ in data.columns):
        raise ValueError("Incomplete OOF/group match")
    y, lo, hi, raw = (data[c].to_numpy() for c in ["y_true", "conf_low", "conf_high", "prediction"])
    folds = data["fold"].to_numpy()
    unique_groups, inverse = np.unique(data["group"].to_numpy(), return_inverse=True)
    group_rows = [np.flatnonzero(inverse == i) for i in range(len(unique_groups))]
    rng = np.random.default_rng(cfg["seed"])
    bootstrap_indices = [
        np.concatenate(
            [group_rows[i] for i in rng.integers(0, len(unique_groups), len(unique_groups))]
        )
        for _ in range(cfg["bootstrap_repeats"])
    ]
    fold_rows = []
    bootstrap_rows = []
    summary = {
        "n": data.height,
        "n_groups": len(unique_groups),
        "bootstrap_repeats": cfg["bootstrap_repeats"],
        "candidates": {},
    }
    for candidate in cfg["candidates"]:
        name = candidate["name"]
        adjusted = raw.mean() + candidate["shift"] + candidate["scale"] * (raw - raw.mean())
        for fold in np.unique(folds):
            mask = folds == fold
            base_score = st_rae(lo[mask], hi[mask], y[mask], raw[mask])
            candidate_score = st_rae(lo[mask], hi[mask], y[mask], adjusted[mask])
            fold_rows.append(
                {
                    "candidate": name,
                    "fold": int(fold),
                    "n": int(mask.sum()),
                    "base_ST_RAE": base_score,
                    "candidate_ST_RAE": candidate_score,
                    "delta_ST_RAE": candidate_score - base_score,
                }
            )
        changes = []
        for repeat, idx in enumerate(bootstrap_indices):
            delta = st_rae(lo[idx], hi[idx], y[idx], adjusted[idx]) - st_rae(
                lo[idx], hi[idx], y[idx], raw[idx]
            )
            bootstrap_rows.append({"candidate": name, "repeat": repeat, "delta_ST_RAE": delta})
            changes.append(delta)
        summary["candidates"][name] = {
            "shift": candidate["shift"],
            "scale": candidate["scale"],
            "overall_delta_ST_RAE": st_rae(lo, hi, y, adjusted) - st_rae(lo, hi, y, raw),
            "folds_improved": sum(
                row["delta_ST_RAE"] < 0 for row in fold_rows if row["candidate"] == name
            ),
            "bootstrap_delta_95pct": np.quantile(changes, [0.025, 0.975]).tolist(),
            "bootstrap_fraction_improved": float(np.mean(np.array(changes) < 0)),
        }
    pl.DataFrame(fold_rows).write_csv(out / "fold-metrics.csv")
    pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
