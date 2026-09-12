"""Combine fixed calibration screen metrics without fitting endpoint selections."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


a = (
    pl.read_csv(path(0))
    .filter(pl.col("fold") == -1)
    .drop("fold")
    .with_columns(pl.lit("affine").alias("family"))
)
b = (
    pl.read_csv(path(1))
    .rename({"base_method": "base"})
    .with_columns(pl.lit("nonlinear").alias("family"))
)
cols = ["base", "calibration", "family", "endpoint", "n", "ST_RAE", "MAE", "Spearman"]
s = pl.concat([a.select(cols), b.select(cols)])
reference = a.filter(pl.col("calibration") == "identity").select(
    "base", "endpoint", pl.col("ST_RAE").alias("baseline")
)
s = s.join(reference, on=["base", "endpoint"], validate="m:1").with_columns(
    (pl.col("ST_RAE") - pl.col("baseline")).alias("delta_ST_RAE")
)
s.sort("base", "endpoint", "ST_RAE").write_csv(out / "all-metrics.csv")
s.group_by("base", "family", "calibration").agg(
    pl.col("ST_RAE").mean().alias("macro"), pl.col("delta_ST_RAE").mean().alias("delta")
).sort("base", "macro").write_csv(out / "all-macro.csv")
s.filter(pl.col("endpoint").str.contains("2D6")).sort("base", "ST_RAE").write_csv(
    out / "d6-comparison.csv"
)
(out / "summary.json").write_text(
    json.dumps(
        dict(
            nominal_methods=31,
            unique_methods=30,
            duplicate="identity is present in both screens",
            readout="Fixed method comparisons only; no endpoint-wise winner assembled",
        ),
        indent=2,
    )
)
