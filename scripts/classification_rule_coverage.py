"""Describe supplied-label coverage without replacing official labels."""

import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
table = pl.read_parquet(next(inp.glob("*-simple-rule-comparison.parquet")))
table = table.filter(pl.col("raw_label").is_not_null()).with_columns(
    pl.col("direct").is_finite().fill_null(False).alias("direct_available"),
    pl.col("tdi").is_finite().fill_null(False).alias("tdi_available"),
    pl.col("simple_label").is_not_null().alias("rule_defined"),
)
summary = (
    table.group_by("endpoint", "raw_label", "direct_available", "tdi_available", "rule_defined")
    .len()
    .sort("endpoint", "raw_label", "direct_available", "tdi_available")
)
summary.write_csv(out / "coverage.csv")
table.filter(~pl.col("rule_defined")).write_csv(out / "unresolved-supplied-labels.csv")
print(summary)
