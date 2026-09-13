"""Compare two validated regression submission vectors without target labels."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

old = pl.read_parquet(inp / cfg["previous"]).sort("Molecule_Name")
new = pl.read_parquet(inp / cfg["current"]).sort("Molecule_Name")
assert old.select("SMILES", "Molecule_Name").equals(new.select("SMILES", "Molecule_Name"))
endpoints = [c for c in old.columns if c.startswith("CYP")]
assert endpoints == [c for c in new.columns if c.startswith("CYP")]

rows = []
for endpoint in endpoints:
    before = old[endpoint].to_numpy()
    after = new[endpoint].to_numpy()
    difference = after - before
    rows.append(
        {
            "endpoint": endpoint,
            "previous_mean": float(before.mean()),
            "current_mean": float(after.mean()),
            "mean_shift": float(difference.mean()),
            "previous_std": float(before.std()),
            "current_std": float(after.std()),
            "pearson": float(np.corrcoef(before, after)[0, 1]),
            "mean_absolute_difference": float(np.mean(np.abs(difference))),
            "root_mean_squared_difference": float(np.sqrt(np.mean(difference**2))),
            "max_absolute_difference": float(np.max(np.abs(difference))),
            "fraction_abs_difference_gt_0_1": float(np.mean(np.abs(difference) > 0.1)),
            "fraction_abs_difference_gt_0_25": float(np.mean(np.abs(difference) > 0.25)),
            "fraction_abs_difference_gt_0_5": float(np.mean(np.abs(difference) > 0.5)),
        }
    )

summary = pl.DataFrame(rows)
summary.write_csv(out / "endpoint-comparison.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "previous_submission": cfg["previous_submission"],
            "current_candidate": cfg["current_candidate"],
            "rows": old.height,
            "endpoints": rows,
        },
        indent=2,
    )
)
