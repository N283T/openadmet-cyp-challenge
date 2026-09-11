"""Prepare five source-specific -Fit_LogAC50 targets for scratch pretraining."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(n):
    return pl.read_parquet(next(inp.glob("*-" + n + ".parquet")))


r = read("pubchem-curated-records").with_columns(
    (-pl.col("Fit_LogAC50").cast(pl.Float64, strict=False)).alias("pac50")
)
accepted = r.filter(
    pl.col("connectivity").is_not_null()
    & ~pl.col("challenge_overlap")
    & pl.col("pac50").is_finite()
)
assert set(accepted["connectivity"]).isdisjoint(set(read("challenge-connectivity")["connectivity"]))
accepted.write_parquet(out / "selected-records.parquet")
targets = accepted.group_by("compound_id", "isoform").agg(
    pl.col("pac50").mean().alias("value"), pl.col("PUBCHEM_SID").n_unique().alias("n_sids")
)
targets.write_parquet(out / "targets.parquet")
compounds = (
    accepted.select("compound_id", "standardized_smiles", "connectivity")
    .unique()
    .sort("compound_id")
)
assert compounds.height == compounds["compound_id"].n_unique()
compounds.write_parquet(out / "external-compounds.parquet")
coverage = r.group_by("isoform", "PUBCHEM_ACTIVITY_OUTCOME").agg(
    pl.len().alias("raw_rows"),
    pl.col("pac50").is_finite().sum().alias("finite_pac50"),
    pl.col("pac50").min().alias("minimum"),
    pl.col("pac50").max().alias("maximum"),
)
coverage.write_csv(out / "raw-coverage.csv")
counts = targets.group_by("isoform").agg(
    pl.len().alias("n"),
    pl.col("value").min().alias("min"),
    pl.col("value").max().alias("max"),
    pl.col("value").mean().alias("mean"),
)
counts.write_csv(out / "target-counts.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            raw_rows=r.height,
            selected_rows=accepted.height,
            compounds=compounds.height,
            labels=targets.height,
            counts=counts.to_dicts(),
            target="negative Fit_LogAC50; original definition log10 molar AC50",
            filter="finite recorded fits, all original outcomes; no invented inactive potency",
            aggregation="mean per standardized compound and CYP",
            exclusion="all challenge connectivity blocks, including test",
        ),
        indent=2,
    )
)
print(counts)
