"""Inspect declared selection-diagnostic inputs."""

import os
from pathlib import Path

import polars as pl

from cyp_challenge.storage import write_json

rows = {}
for p in Path(os.environ["CYP_INPUT_DIR"]).glob("*.parquet"):
    f = pl.read_parquet(p)
    rows[p.name] = {
        "schema": {k: str(v) for k, v in f.schema.items()},
        "head": f.head(2).to_dicts(),
    }
write_json(Path(os.environ["CYP_OUTPUT_DIR"]) / "schema.json", rows)
print(rows)
