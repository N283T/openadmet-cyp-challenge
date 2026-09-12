"""Merge source-separated external assays and mask contradictory binary supervision."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(
        next(p for p in inp.iterdir() if p.name.partition("-")[2] == name + ".parquet")
    )


c, o, t, tc = (
    read("external-compounds"),
    read("observations"),
    read("nvs-targets"),
    read("nvs-compounds"),
)
t.write_parquet(out / "nvs-source-targets.parquet")
raw = read("nvs-curated-records")
raw.write_parquet(out / "nvs-source-records.parquet")
raw.select("aeid", "aenm", "tested_conc_unit", "resp_unit").unique().write_parquet(
    out / "nvs-assay-metadata.parquet"
)
t = t.join(tc.select("compound_id", "connectivity"), on="compound_id", validate="m:1")
new = t.select(
    "compound_id",
    "isoform",
    "connectivity",
    pl.lit(False).alias("challenge_overlap"),
    (
        pl.lit("nvs:")
        + pl.col("aeid").cast(pl.String)
        + ":"
        + pl.col("compound_id")
        + ":"
        + pl.col("label_kind")
    ).alias("observation_id"),
    "source",
    (pl.lit("nvs:") + pl.col("aeid").cast(pl.String) + ":" + pl.col("isoform") + ":AC50").alias(
        "assay_id"
    ),
    pl.lit("AC50").alias("measurement_type"),
    pl.when(pl.col("label_kind") == "pac50").then(pl.col("value")).alias("value"),
    pl.when(pl.col("label_kind") == "pac50").then(pl.lit("=")).alias("relation"),
    pl.when(pl.col("label_kind") == "pac50").then(pl.col("value")).alias("lower_bound"),
    pl.when(pl.col("label_kind") == "pac50").then(pl.col("value")).alias("upper_bound"),
    pl.when(pl.col("label_kind") == "binary").then(pl.col("value")).alias("binary_label"),
    pl.when(pl.col("label_kind") == "pac50")
    .then(pl.lit("Active"))
    .when(pl.col("value") == 1)
    .then(pl.lit("Active"))
    .when(pl.col("value") == 0)
    .then(pl.lit("Inactive"))
    .otherwise(pl.lit("Conflicting"))
    .alias("hitcall"),
    pl.lit(None, dtype=pl.Float64).alias("max_response"),
    pl.lit(None, dtype=pl.Float64).alias("fit_r2"),
    pl.lit(None, dtype=pl.Float64).alias("curve_class"),
    pl.when(pl.col("label_kind") == "pac50")
    .then(pl.lit("hitc=1; representative gsid_rep=1; uM; finite unflagged; pAC50=6-modl_ga"))
    .otherwise(pl.lit("hitc=0/1; representative gsid_rep=1; uM; duplicate means"))
    .alias("qc"),
    pl.lit(False).alias("preincubation_mentioned"),
)
combined = (
    pl.concat([o, new], how="diagonal_relaxed")
    .with_columns(
        pl.coalesce("original_binary_label", "binary_label").alias("original_binary_label")
    )
    .drop("binary_conflict")
)
conflicts = combined.group_by("compound_id", "assay_id").agg(
    (
        (pl.col("original_binary_label").drop_nulls().n_unique() > 1)
        | (~pl.col("original_binary_label").drop_nulls().is_in([0.0, 1.0])).any()
    ).alias("binary_conflict")
)
combined = combined.join(conflicts, on=["compound_id", "assay_id"], validate="m:1").with_columns(
    pl.when(~pl.col("binary_conflict")).then(pl.col("binary_label")).alias("binary_label")
)
compounds = pl.concat([c, tc]).unique().sort("compound_id")
assert compounds["compound_id"].n_unique() == compounds.height
assert combined["observation_id"].n_unique() == combined.height
assert combined.join(compounds, on="compound_id", how="anti").height == 0
assert combined["binary_label"].drop_nulls().is_in([0.0, 1.0]).all()
assert set(compounds["connectivity"]).isdisjoint(read("challenge-connectivity")["connectivity"])
assert not combined["challenge_overlap"].any()
combined.sort("observation_id").write_parquet(out / "observations.parquet")
compounds.write_parquet(out / "external-compounds.parquet")
conflicts.filter("binary_conflict").write_parquet(out / "binary-conflicts.parquet")
inventory = combined.group_by("source", "isoform", "measurement_type").agg(
    pl.len().alias("observations"),
    pl.col("compound_id").n_unique().alias("compounds"),
    pl.col("assay_id").n_unique().alias("assays"),
    pl.col("value").is_not_null().sum().alias("continuous"),
    pl.col("binary_label").is_not_null().sum().alias("binary"),
)
inventory.sort("source", "isoform").write_csv(out / "inventory.csv")
summary = dict(
    compounds=compounds.height,
    new_compounds=compounds.height - c.height,
    observations=combined.height,
    continuous=combined["value"].is_not_null().sum(),
    binary=combined["binary_label"].is_not_null().sum(),
    masked_conflicting_compound_assays=conflicts.filter("binary_conflict").height,
    binary_policy=(
        "Contradictory percompound/assay binary labels and fractional duplicate means masked; "
        "original_binary_label retained"
    ),
    source_policy=(
        "Different source assays retain separate heads and records; "
        "exact compound identity deduplicated only"
    ),
    nvs_provenance=(
        "Original NVS observations, concentration units and QC remain in source run "
        "20260912T052128343558Z-128af2eb; aggregated targets also copied"
    ),
    inventory=inventory.to_dicts(),
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
