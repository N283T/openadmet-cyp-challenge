"""Materialize anchor-retained memberships without splitting connectivity blocks."""

import os
from pathlib import Path

import polars as pl
from rdkit import Chem

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(n):
    return pl.read_parquet(next(inp.glob("*-" + n + ".parquet")))


a = read("assignments")
c = read("compounds")
labels = read("measurements")
s = pl.read_csv(next(inp.glob("*-hit-seeds.csv")))
c = c.with_columns(
    pl.Series(
        "connectivity",
        [Chem.MolToInchiKey(Chem.MolFromSmiles(x)).split("-")[0] for x in c["standardized_smiles"]],
    )
)
anchor_keys = set(c.filter(pl.col("compound_id").is_in(s["compound_id"].to_list()))["connectivity"])
a = a.join(c.select("compound_id", "connectivity"), on="compound_id", validate="1:1").with_columns(
    pl.col("connectivity").is_in(list(anchor_keys)).alias("is_anchor")
)
a = a.with_columns(((pl.col("hit_origins") != "") & ~pl.col("is_anchor")).alias("is_neighbor"))
a.write_parquet(out / "parent-assignments.parquet")
rows = []
for fold in range(5):
    f = a.with_columns(
        pl.when((pl.col("fold") == fold) & ~pl.col("is_anchor"))
        .then(pl.lit("validation"))
        .otherwise(pl.lit("train"))
        .alias("role"),
        pl.lit(fold).alias("outer_fold"),
    )
    assert not f.filter(pl.col("is_anchor") & (pl.col("role") != "train")).height
    assert f.group_by("connectivity").agg(pl.col("role").n_unique())["role"].max() == 1
    rows.append(f)
m = pl.concat(rows)
m.write_parquet(out / "membership.parquet")
f = labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")).join(
    a, on="compound_id", validate="m:1"
)
f.group_by("endpoint").agg(
    pl.len().alias("total"),
    pl.col("is_anchor").sum().alias("anchors"),
    pl.col("is_neighbor").sum().alias("neighbors"),
).write_csv(out / "coverage.csv")
write_json(
    out / "summary.json",
    dict(
        anchor_compounds=a["is_anchor"].sum(),
        neighbor_compounds=a["is_neighbor"].sum(),
        anchor_always_train=True,
        connectivity_integrity=True,
        aux_policy="Auxiliary excludes full fold pairs; direct readouts retain anchors.",
    ),
)
