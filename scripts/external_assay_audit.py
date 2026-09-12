"""Audit usable labels, assay sizes, overlapping labels and identity exclusions."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(
        next(p for p in inp.iterdir() if p.name.partition("-")[2] == name + ".parquet")
    )


o, c, a = read("observations"), read("external-compounds"), read("assay-inventory")
assert o["observation_id"].n_unique() == o.height
assert o.join(c, on="compound_id", how="anti").height == 0
assert not o["challenge_overlap"].any()
assert (
    not o.filter((pl.col("source") == "pubchem1851") & (pl.col("hitcall") == "Inactive"))["value"]
    .is_not_null()
    .any()
)
assert not o.filter(pl.col("relation") != "=")["value"].is_not_null().any()
assert set(c["connectivity"]).isdisjoint(read("challenge-connectivity")["connectivity"])
conflicts = o.group_by("compound_id", "assay_id").agg(
    pl.len().alias("records"),
    pl.col("binary_label").drop_nulls().n_unique().alias("distinct_binary"),
    pl.col("value").max().alias("max_value"),
    pl.col("value").min().alias("min_value"),
)
conflicts.filter(pl.col("records") > 1).write_parquet(out / "replicate-audit.parquet")
g = (
    o.filter(pl.col("value").is_not_null())
    .group_by("assay_id")
    .agg(
        pl.col("compound_id").n_unique().alias("n_compounds"),
        pl.col("value").n_unique().alias("n_values"),
    )
)
summary = dict(
    compounds=c.height,
    observations=o.height,
    continuous=o["value"].is_not_null().sum(),
    binary=o["binary_label"].is_not_null().sum(),
    censored=o["relation"].is_in(["<", ">", "<=", ">="]).sum(),
    ranking_assays_2=g.filter((pl.col("n_compounds") >= 2) & (pl.col("n_values") >= 2)).height,
    ranking_assays_5=g.filter((pl.col("n_compounds") >= 5) & (pl.col("n_values") >= 2)).height,
    classification_assays_5=a.filter(pl.col("binary") >= 5).height,
    replicate_keys=conflicts.filter(pl.col("records") > 1).height,
    conflicting_binary_keys=conflicts.filter(pl.col("distinct_binary") > 1).height,
)
qc = (
    o.filter(pl.col("source") == "pubchem1851")
    .group_by("isoform", "hitcall", "curve_class")
    .agg(
        pl.len().alias("n"),
        pl.col("value").is_not_null().sum().alias("finite_pac50"),
        pl.col("fit_r2").median().alias("median_fit_r2"),
    )
)
qc.write_csv(out / "pubchem-curve-quality.csv")
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
