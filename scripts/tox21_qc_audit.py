"""Audit Tox21 phenotype conflicts and replicate usability before model use."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
r = pl.read_parquet(next(inp.glob("*-tox21-curated-replicates.parquet")))
t = pl.read_parquet(next(inp.glob("*-tox21-targets.parquet")))
r.group_by("isoform", "PUBCHEM_ACTIVITY_OUTCOME", "Phenotype").agg(
    pl.len().alias("rows"), pl.col("pac50").is_finite().sum().alias("finite_pac50")
).sort("isoform", "PUBCHEM_ACTIVITY_OUTCOME", "Phenotype").write_csv(out / "phenotype-audit.csv")
b = t.filter(pl.col("label_kind") == "binary")
b.group_by("isoform").agg(
    pl.len().alias("n"),
    ((pl.col("value") > 0) & (pl.col("value") < 1)).sum().alias("conflicting_compounds"),
).write_csv(out / "binary-conflicts.csv")
for col in ["Compound QC", "Analysis Comment", "Curve_Description"]:
    r.group_by("isoform", col).len().sort("isoform", col).write_csv(
        out / (col.replace(" ", "-") + ".csv")
    )
summary = (
    r.filter(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active")
    .group_by("isoform")
    .agg(
        (pl.col("Phenotype") == "Activator").sum().alias("active_activator_replicates"),
        (pl.col("Phenotype") == "Inhibitor").sum().alias("active_inhibitor_replicates"),
    )
)
(out / "summary.json").write_text(json.dumps(summary.to_dicts(), indent=2))
print(summary, flush=True)
