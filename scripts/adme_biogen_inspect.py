"""Inspect author CSV field definitions and original censor tokens."""

import json
import os
from pathlib import Path

import polars as pl

p = next(Path(os.environ["CYP_INPUT_DIR"]).glob("*-ADME_public_set_3521.csv"))
d = pl.read_csv(p, infer_schema=False)
r = dict(
    rows=d.height,
    columns=d.columns,
    examples=d.head(3).to_dicts(),
    missing={c: d[c].null_count() for c in d.columns},
    censored={
        c: d.filter(pl.col(c).str.contains(r"^[<>~=]"))
        .select(c)
        .unique()
        .to_series()
        .to_list()[:20]
        for c in d.columns
    },
)
Path(os.environ["CYP_OUTPUT_DIR"], "inspection.json").write_text(json.dumps(r, indent=2))
print(json.dumps(r, indent=2))
