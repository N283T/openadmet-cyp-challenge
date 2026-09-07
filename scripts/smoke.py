"""Deterministic infrastructure probe; run via `cyp run`, never as an experiment template."""

import os
from pathlib import Path

import polars as pl

source = Path(os.environ["CYP_INPUT_DIR"])
output = Path(os.environ["CYP_OUTPUT_DIR"])
rows = []
for path in sorted(source.glob("*.csv")):
    table = pl.read_csv(path)
    rows.append({"file": path.name, "rows": table.height, "columns": table.width})
pl.DataFrame(rows).write_json(output / "summary.json")
print("Persisted input summary:", rows)
