"""Audit observable label availability, primary response and chemical coverage."""

import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


comp, aliases, measure, aux = [read(n) for n in ["compounds", "aliases", "measurements", "targets"]]
testids = set(aliases.filter(pl.col("role") == "test")["compound_id"])
comp = comp.sort("compound_id")
ids = comp["compound_id"].to_list()
mols = [Chem.MolFromSmiles(s) for s in comp["standardized_smiles"]]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(m) for m in mols]
lookup = {c: i for i, c in enumerate(ids)}
desc = comp.select("compound_id").with_columns(
    pl.Series("MW", [Descriptors.MolWt(m) for m in mols]),
    pl.Series("logP", [Descriptors.MolLogP(m) for m in mols]),
    pl.Series("TPSA", [Descriptors.TPSA(m) for m in mols]),
)
rows, stats, corr, bins, neigh, coverage = [], [], [], [], [], []
for enzyme in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
    ep = enzyme + "_pIC50_direct_inhibition"
    labels = measure.filter(pl.col("endpoint") == ep).select(
        "compound_id", pl.col("value").alias("pic50")
    )
    primary = aux.filter(pl.col("head").str.starts_with(enzyme + "@")).select(
        "compound_id", pl.col("value").alias("log2fc")
    )
    assert primary["compound_id"].n_unique() == primary.height
    f = (
        comp.select("compound_id")
        .join(labels, on="compound_id", how="left", validate="1:1")
        .join(primary, on="compound_id", how="left", validate="1:1")
        .join(desc, on="compound_id", validate="1:1")
    )
    f = f.with_columns(
        pl.col("compound_id").is_in(list(testids)).alias("is_test"),
        pl.col("pic50").is_not_null().alias("has_label"),
        pl.lit(enzyme).alias("enzyme"),
    )
    f = f.with_columns(
        pl.when(pl.col("is_test"))
        .then(pl.lit("test"))
        .when(pl.col("has_label"))
        .then(pl.lit("labeled"))
        .otherwise(pl.lit("unlabeled"))
        .alias("population")
    )
    rows.append(f)
    for group, g in f.group_by("population"):
        for variable in ["log2fc", "pic50", "MW", "logP", "TPSA"]:
            x = g[variable].drop_nulls().to_numpy()
            if len(x):
                stats.append(
                    dict(
                        enzyme=enzyme,
                        population=group[0],
                        variable=variable,
                        n=len(x),
                        mean=float(x.mean()),
                        median=float(np.median(x)),
                        q10=float(np.quantile(x, 0.1)),
                        q90=float(np.quantile(x, 0.9)),
                    )
                )
        coverage.append(
            dict(enzyme=enzyme, population=group[0], n=g.height, primary_n=g["log2fc"].count())
        )
    p = f.filter(~pl.col("is_test") & pl.col("log2fc").is_not_null())
    both = p.filter(pl.col("has_label"))
    corr.append(
        dict(
            enzyme=enzyme,
            primary_n=p.height,
            primary_labeled_n=both.height,
            availability_auc=float(roc_auc_score(p["has_label"], -p["log2fc"])),
            response_pic50_spearman=float(spearmanr(both["log2fc"], both["pic50"]).statistic),
        )
    )
    for low, high in [(-100, -2), (-2, -1), (-1, -0.5), (-0.5, 0), (0, 100)]:
        g = p.filter((pl.col("log2fc") >= low) & (pl.col("log2fc") < high))
        bins.append(
            dict(
                enzyme=enzyme,
                low=low,
                high=high,
                n=g.height,
                labeled_n=int(g["has_label"].sum()),
                label_fraction=float(g["has_label"].mean()) if g.height else None,
            )
        )
    pool = [lookup[c] for c in labels["compound_id"]]
    for group in ["test", "unlabeled"]:
        for c in f.filter(pl.col("population") == group)["compound_id"]:
            sims = np.asarray(
                DataStructs.BulkTanimotoSimilarity(fps[lookup[c]], [fps[i] for i in pool])
            )
            j = int(sims.argmax())
            neigh.append(
                dict(
                    enzyme=enzyme,
                    population=group,
                    compound_id=c,
                    nearest_labeled_id=ids[pool[j]],
                    tanimoto=float(sims[j]),
                )
            )
pl.concat(rows).write_parquet(out / "availability.parquet")
for name, data in [
    ("distributions", stats),
    ("coverage", coverage),
    ("associations", corr),
    ("primary-bins", bins),
]:
    pl.DataFrame(data).write_csv(out / (name + ".csv"))
n = pl.DataFrame(neigh)
n.write_parquet(out / "nearest-labeled.parquet")
n.group_by("enzyme", "population").agg(
    pl.col("tanimoto").median().alias("median"),
    pl.col("tanimoto").quantile(0.1).alias("q10"),
    pl.col("tanimoto").quantile(0.9).alias("q90"),
).write_csv(out / "similarity-summary.csv")
write_json(
    out / "summary.json",
    dict(
        associations=corr,
        coverage=coverage,
        limits=(
            "Public label presence is not a selection/failed-fit flag. "
            "Test responses are unavailable; no causal or hidden-label inference."
        ),
    ),
)
print(pl.DataFrame(corr))
