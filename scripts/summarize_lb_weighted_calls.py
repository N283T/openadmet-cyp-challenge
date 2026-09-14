"""Summarize deployment call counts for the LB-prior-weighted screen."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
deployment = pl.read_parquet(inp / cfg["deployment"])
summary = (
    deployment.group_by("endpoint", "method")
    .agg(pl.col("binary").sum().alias("positives"), pl.col("threshold").first())
    .sort("endpoint", "positives")
)
summary.write_csv(out / "call-counts.csv")
print(summary.write_csv(), flush=True)
