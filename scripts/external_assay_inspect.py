"""Inspect immutable source schemas for assay-aware external preparation."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
result = {}
for p in inp.iterdir():
    if p.suffix != ".parquet":
        continue
    d = pl.read_parquet(p)
    result[p.name] = {"shape": d.shape, "schema": str(d.schema), "sample": d.head(2).to_dicts()}
    if "standard_relation" in d.columns:
        result[p.name]["relations"] = (
            d.group_by("standard_relation", "standard_units", "standard_type").len().to_dicts()
        )
Path(os.environ["CYP_OUTPUT_DIR"], "inspection.json").write_text(
    json.dumps(result, default=str, indent=2)
)
print(json.dumps(result, default=str, indent=2))
