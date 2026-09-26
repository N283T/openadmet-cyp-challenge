"""Paired comparison of the original squared and robust D6-bound blends."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def score(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    return float(
        (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
        / (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold", "split_id"]
    old = (
        pl.read_parquet(inp / "squared-blend.parquet")
        .filter(
            (pl.col("constraint_weight") == cfg["squared_constraint_weight"])
            & (pl.col("placement") == "fold_moment_match")
            & (pl.col("blend_weight") == 0.5)
        )
        .select(
            *keys, "y_true", "conf_low", "conf_high", pl.col("prediction").alias("old_prediction")
        )
    )
    new = (
        pl.read_parquet(inp / "robust-blend.parquet")
        .filter(
            (pl.col("robust_delta") == cfg["robust_delta"])
            & (pl.col("constraint_weight") == cfg["robust_constraint_weight"])
            & (pl.col("placement") == "fold_moment_match")
            & (pl.col("blend_weight") == 0.5)
        )
        .select(
            *keys, "y_true", "conf_low", "conf_high", pl.col("prediction").alias("new_prediction")
        )
    )
    if old.height != 2914 or new.height != 2914:
        raise ValueError("Incomplete candidate OOF")
    paired = old.join(new, on=keys, validate="1:1", suffix="_new")
    if paired.height != old.height:
        raise ValueError("Unmatched candidates")
    for key in ["y_true", "conf_low", "conf_high"]:
        if not np.array_equal(paired[key].to_numpy(), paired[f"{key}_new"].to_numpy()):
            raise ValueError(f"Mismatched {key}")
    summary = []
    for split in ["parent", "scaffold"]:
        part = paired.filter(pl.col("split_id") == split).join(
            pl.read_parquet(inp / f"{split}-assignments.parquet").select(
                "compound_id", "fold", "group"
            ),
            on=["compound_id", "fold"],
            validate="1:1",
        )
        y, lo, hi, p0, p1 = [
            part[col].to_numpy()
            for col in ["y_true", "conf_low", "conf_high", "old_prediction", "new_prediction"]
        ]
        unique, inverse = np.unique(part["group"].to_numpy(), return_inverse=True)
        members = [np.flatnonzero(inverse == j) for j in range(len(unique))]
        rng = np.random.default_rng(cfg["seed"] + (split == "scaffold"))
        draws = np.empty(cfg["bootstrap_draws"])
        for j in range(len(draws)):
            ids = np.concatenate([members[k] for k in rng.integers(len(unique), size=len(unique))])
            draws[j] = score(y[ids], lo[ids], hi[ids], p1[ids]) - score(
                y[ids], lo[ids], hi[ids], p0[ids]
            )
        summary.append(
            {
                "split_id": split,
                "n": part.height,
                "squared_ST_RAE": score(y, lo, hi, p0),
                "robust_ST_RAE": score(y, lo, hi, p1),
                "robust_minus_squared": score(y, lo, hi, p1) - score(y, lo, hi, p0),
                "q025": float(np.quantile(draws, 0.025)),
                "q975": float(np.quantile(draws, 0.975)),
                "p_robust_better": float(np.mean(draws < 0)),
            }
        )
    pl.DataFrame(summary).write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")


if __name__ == "__main__":
    main()
