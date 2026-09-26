"""Grouped paired bootstrap for a fixed CYP1A2 atomic ensemble contrast."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    error = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    mean = y.mean()
    baseline = (np.maximum(lo - mean, 0) + np.maximum(mean - hi, 0)).sum()
    return float(error / baseline)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    oof = pl.read_parquet(inp / "component-oof.parquet").filter(
        pl.col("endpoint") == cfg["endpoint"]
    )
    rows = []
    for split in ["parent", "scaffold"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "group", "is_anchor"
        )
        frame = oof.filter(pl.col("split_id") == split).join(
            assignments, on=["compound_id", "fold"], validate="m:1"
        )
        if frame.height != oof.filter(pl.col("split_id") == split).height:
            raise ValueError(f"Incomplete {split} groups")
        if frame["is_anchor"].any():
            raise ValueError("Anchor appeared in validation OOF")
        wide = frame.pivot(
            on="method", index=["compound_id", "endpoint", "fold", "group",
                                "y_true", "conf_low", "conf_high"], values="prediction"
        ).sort("compound_id")
        if wide.height != 1372:
            raise ValueError(f"Incomplete {split} pivot")
        y = wide["y_true"].to_numpy()
        lo = wide["conf_low"].to_numpy()
        hi = wide["conf_high"].to_numpy()
        target = wide[cfg["target_recipe"]].to_numpy()
        groups = wide["group"].to_numpy()
        unique_groups = np.unique(groups)
        group_indices = [np.flatnonzero(groups == group) for group in unique_groups]
        rng = np.random.default_rng(cfg["seed"] + (0 if split == "parent" else 1))
        for reference in cfg["references"]:
            control = wide[reference].to_numpy()
            observed = st_rae(y, lo, hi, target) - st_rae(y, lo, hi, control)
            deltas = np.empty(cfg["replicates"], dtype=float)
            for replicate in range(cfg["replicates"]):
                sampled = rng.integers(len(group_indices), size=len(group_indices))
                index = np.concatenate([group_indices[i] for i in sampled])
                deltas[replicate] = (
                    st_rae(y[index], lo[index], hi[index], target[index])
                    - st_rae(y[index], lo[index], hi[index], control[index])
                )
            rows.append({
                "split_id": split, "reference": reference,
                "target": cfg["target_recipe"], "n": len(y),
                "groups": len(unique_groups), "observed_delta": observed,
                "lower_95": float(np.quantile(deltas, 0.025)),
                "upper_95": float(np.quantile(deltas, 0.975)),
                "fraction_better": float((deltas < 0).mean()),
            })
    pl.DataFrame(rows).sort("split_id", "reference").write_csv(out / "bootstrap.csv")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg, "comparisons": len(rows),
        "qualification": "Conditional group bootstrap; does not adjust for recipe selection",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
