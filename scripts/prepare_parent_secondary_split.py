"""Build an anchor-retained secondary split from an immutable base assignment."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.storage import write_json

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

base = pl.read_parquet(inp / "base-assignments.parquet").filter(
    pl.col("strategy") == cfg["strategy"]
)
parent = pl.read_parquet(inp / "primary-parent-assignments.parquet")
metadata_columns = [
    "compound_id",
    "hit_origins",
    "connectivity",
    "is_anchor",
    "is_neighbor",
]

if base["compound_id"].n_unique() != base.height:
    raise ValueError("Secondary assignment has duplicate compound IDs")
if parent["compound_id"].n_unique() != parent.height:
    raise ValueError("Primary parent assignment has duplicate compound IDs")

assignments = base.join(parent.select(metadata_columns), on="compound_id", validate="1:1").sort(
    "compound_id"
)
if assignments.height != parent.height:
    raise ValueError("Secondary assignment does not cover the primary parent population")
if assignments.group_by("connectivity").agg(pl.col("fold").n_unique().alias("n"))["n"].max() != 1:
    raise ValueError("A connectivity block crosses secondary folds")
if sorted(assignments["fold"].unique().to_list()) != list(range(cfg["folds"])):
    raise ValueError("Secondary assignment does not contain the expected folds")

assignments.write_parquet(out / "parent-assignments.parquet")
memberships = []
for fold in range(cfg["folds"]):
    membership = assignments.with_columns(
        pl.when((pl.col("fold") == fold) & ~pl.col("is_anchor"))
        .then(pl.lit("validation"))
        .otherwise(pl.lit("train"))
        .alias("role"),
        pl.lit(fold).alias("outer_fold"),
    )
    if membership.filter(pl.col("is_anchor") & (pl.col("role") != "train")).height:
        raise ValueError("An anchor was withheld from training")
    memberships.append(membership)
pl.concat(memberships).write_parquet(out / "membership.parquet")

fold_counts = (
    assignments.group_by("fold")
    .agg(
        pl.len().alias("all_compounds"),
        pl.col("is_anchor").sum().alias("anchors"),
        (~pl.col("is_anchor")).sum().alias("validation_compounds"),
    )
    .sort("fold")
)
fold_counts.write_csv(out / "fold-counts.csv")
write_json(
    out / "summary.json",
    {
        "strategy": cfg["strategy"],
        "role": "secondary_ordering_check",
        "selection_policy": cfg["selection_policy"],
        "compounds": assignments.height,
        "anchors": int(assignments["is_anchor"].sum()),
        "connectivity_integrity": True,
        "anchor_always_train": True,
        "fold_counts": fold_counts.to_dicts(),
    },
)
