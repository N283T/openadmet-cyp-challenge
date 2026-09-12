"""Verify the serialized assay exclusion regex against explicit TDI cases."""

import json
import os
import re
from pathlib import Path

import polars as pl

cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
pattern = cfg["tdi_regex"]
cases = [
    "TDI",
    "Assay assessing TDI inhibition",
    "time-dependent inhibition",
    "irreversible inhibition",
    "ordinary inhibition",
]
result = {s: bool(re.search(pattern, s)) for s in cases}
p = next(Path(os.environ["CYP_INPUT_DIR"]).glob("*-chembl-raw.parquet"))
r = pl.read_parquet(p)
result["matched_source_rows"] = r.filter(pl.col("assay_description").str.contains(pattern)).height
result["explicit_tdi_source_rows"] = r.filter(
    pl.col("assay_description").str.contains(r"(?i)\bTDI\b")
).height
Path(os.environ["CYP_OUTPUT_DIR"], "regex-audit.json").write_text(json.dumps(result, indent=2))
print(result)
assert result["TDI"] and result["Assay assessing TDI inhibition"]
assert not result["ordinary inhibition"]
