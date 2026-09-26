"""Place the completed D6 blend at the historically scored D6 moments."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    old = pl.read_csv(inp / "current.csv").sort("Molecule_Name")
    candidate = pl.read_csv(inp / "candidate.csv").sort("Molecule_Name")
    reference = pl.read_csv(inp / "reference.csv").sort("Molecule_Name")
    if (
        old.height != 750
        or not old.select(reference.columns).equals(reference)
        or not candidate.select("Molecule_Name", "SMILES").equals(
            old.select("Molecule_Name", "SMILES")
        )
    ):
        raise ValueError("Molecule/reference mismatch")
    for col in old.columns[2:]:
        if col != endpoint and not old[col].equals(candidate[col]):
            raise ValueError(f"Changed unrelated endpoint: {col}")
    old_y = old[endpoint].to_numpy()
    initial = candidate[endpoint].to_numpy()
    placed = (initial - initial.mean()) / initial.std(ddof=1) * cfg["target_sample_sd"] + cfg[
        "target_mean"
    ]
    if not np.isfinite(placed).all() or placed.min() < 1:
        raise ValueError("Invalid placed predictions")
    result = candidate.with_columns(pl.Series(endpoint, placed))
    result.write_csv(out / "candidate-sd90.csv")
    result.select("Molecule_Name", "SMILES", endpoint).with_columns(
        pl.lit("blind").alias("split_id"),
        pl.lit(-1).alias("fold"),
        pl.Series("current_prediction", old_y),
        pl.Series("original_candidate", initial),
    ).write_parquet(out / "d6-comparison.parquet")
    summary = {
        "config": cfg,
        "rows": result.height,
        "original_mean": float(initial.mean()),
        "original_sample_sd": float(initial.std(ddof=1)),
        "placed_mean": float(placed.mean()),
        "placed_sample_sd": float(placed.std(ddof=1)),
        "mean_absolute_delta": float(np.abs(placed - initial).mean()),
        "maximum_absolute_delta": float(np.abs(placed - initial).max()),
        "spearman_vs_current": float(spearmanr(placed, old_y).statistic),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
