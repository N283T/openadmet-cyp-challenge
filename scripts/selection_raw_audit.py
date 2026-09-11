"""Validate primary aggregation and visualize observed availability; no target imputation."""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def file(name):
    return next(inp.glob("*-" + name))


a = (
    pl.read_parquet(file("aliases.parquet"))
    .filter(pl.col("role") == "aux")
    .select(
        pl.col("molecule_name").alias("Molecule_Name"),
        pl.col("raw_smiles").alias("SMILES"),
        "compound_id",
    )
    .unique()
)
r = pl.read_csv(file("cyp-challenge-single-concentration-TRAIN.csv")).join(
    a, on=["Molecule_Name", "SMILES"], how="left", validate="m:1"
)
assert r["compound_id"].null_count() == 0
v = pl.read_parquet(file("availability.parquet"))
agg = r.group_by("compound_id", "enzyme").agg(
    pl.col("log2fc_estimate").mean().alias("raw_log2fc"),
    pl.col("log2fc_fdr").min().alias("min_fdr"),
    pl.len().alias("raw_rows"),
)
z = v.join(agg, on=["compound_id", "enzyme"], how="left", validate="1:1")
p = z.filter(pl.col("log2fc").is_not_null())
assert np.allclose(p["log2fc"], p["raw_log2fc"], atol=1e-12)
summary = p.group_by("enzyme", "population").agg(
    pl.len().alias("n"),
    pl.col("log2fc").min().alias("min_log2fc"),
    pl.col("log2fc").max().alias("max_log2fc"),
    (pl.col("min_fdr") < 0.05).mean().alias("fraction_fdr_lt_005"),
    (pl.col("pic50") < 4).sum().alias("labeled_pic50_below4"),
)
summary.write_csv(out / "raw-audit.csv")
fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
for enzyme, ax in zip(["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"], axes.flat, strict=True):
    for pop in ["labeled", "unlabeled"]:
        x = p.filter((pl.col("enzyme") == enzyme) & (pl.col("population") == pop))[
            "log2fc"
        ].to_numpy()
        ax.hist(x, bins=np.linspace(-5, 2, 57), histtype="step", label=f"{pop} (n={len(x)})")
    ax.set(title=enzyme, xlabel="Primary TDI log2 fold change", ylabel="Compound count")
    ax.legend()
fig.savefig(out / "primary-availability.png", dpi=160)
write_json(
    out / "summary.json",
    dict(
        raw_rows=r.height,
        primary_compounds=r["compound_id"].n_unique(),
        aggregation_matches=True,
        scope=(
            "FDR < 0.05 is descriptive, not an official cutoff. "
            "Histograms show [-5,2]; tables retain full range."
        ),
    ),
)
print(summary)
