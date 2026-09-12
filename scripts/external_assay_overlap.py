"""Diagnose external/challenge assay agreement; never fit or apply a CV calibration."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(
        next(p for p in inp.iterdir() if p.name.partition("-")[2] == name + ".parquet")
    )


ext = read("challenge-overlap-observations")
labels = (
    read("measurements")
    .filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
    .select(
        "compound_id",
        pl.col("endpoint").str.split("_").list.first().alias("isoform"),
        pl.col("value").alias("challenge_pic50"),
    )
)
assert labels.select("compound_id", "isoform").n_unique() == labels.height
raw, reps = read("tox21-curated-records"), read("tox21-curated-replicates")
active = raw.filter(
    pl.col("challenge_overlap") & (pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active")
).select("aid", "PUBCHEM_SID")
t = reps.join(active, on=["aid", "PUBCHEM_SID"], validate="m:1").filter(
    (pl.col("Phenotype") == "Inhibitor")
    & pl.col("pac50").is_finite()
    & pl.col("curve_class").is_in([-1.1, -1.2, -2.1, -2.2])
)
t = (
    t.group_by("compound_id", "isoform", "aid", "connectivity")
    .agg(pl.col("pac50").median().alias("value"))
    .with_columns(
        pl.lit("tox21_p450_glo").alias("source"),
        (
            pl.lit("tox21:") + pl.col("aid").cast(pl.String) + ":" + pl.col("isoform") + ":AC50"
        ).alias("assay_id"),
        pl.lit("AC50").alias("measurement_type"),
    )
)
cols = ["compound_id", "connectivity", "isoform", "source", "assay_id", "measurement_type", "value"]
ext = pl.concat([ext.filter(pl.col("value").is_not_null()).select(cols), t.select(cols)])
ext = ext.group_by(cols[:-1]).agg(pl.col("value").median().alias("external_potency"))
paired = ext.join(labels, on=["compound_id", "isoform"], validate="m:1")
paired.write_parquet(out / "exact-identity-pairs.parquet")
conn = labels.join(read("challenge-connectivity"), on="compound_id", validate="m:1")
connectivity_pairs = ext.join(
    conn.rename({"compound_id": "challenge_compound_id"}), on=["connectivity", "isoform"]
).with_columns((pl.col("compound_id") == pl.col("challenge_compound_id")).alias("exact_identity"))
connectivity_pairs.write_parquet(out / "connectivity-pairs.parquet")


def stats(g):
    # Each group already one row per compound/assay. Source summaries may repeat compounds.
    x = g["external_potency"].to_numpy()
    y = g["challenge_pic50"].to_numpy()
    offset = y - x
    r = dict(
        n=g.height,
        compounds=g["compound_id"].n_unique(),
        assays=g["assay_id"].n_unique(),
        mean_offset=float(offset.mean()),
        median_offset=float(np.median(offset)),
        offset_sd=float(np.std(offset)),
        offset_q10=float(np.quantile(offset, 0.1)),
        offset_q90=float(np.quantile(offset, 0.9)),
        mae_identity=float(np.mean(abs(offset))),
        mae_descriptive_median_offset=float(np.mean(abs(offset - np.median(offset)))),
    )
    if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0:
        r.update(
            pearson=float(pearsonr(x, y).statistic),
            spearman=float(spearmanr(x, y).statistic),
            descriptive_slope=float(np.polyfit(x, y, 1)[0]),
        )
    return r


rows = []
for key, g in paired.group_by("source", "isoform", "measurement_type"):
    rows.append(dict(zip(["source", "isoform", "measurement_type"], key, strict=True)) | stats(g))
pl.DataFrame(rows, infer_schema_length=None).sort(
    "source", "isoform", "measurement_type"
).write_csv(out / "source-agreement.csv")
rows = []
for key, g in paired.group_by("source", "isoform", "measurement_type", "assay_id"):
    rows.append(
        dict(zip(["source", "isoform", "measurement_type", "assay_id"], key, strict=True))
        | stats(g)
    )
a = pl.DataFrame(rows, infer_schema_length=None)
a.sort("source", "isoform", "assay_id").write_csv(out / "assay-agreement.csv")
summary = dict(
    exact_identity_pairs=paired.height,
    exact_identity_compounds=paired["compound_id"].n_unique(),
    connectivity_pairs=connectivity_pairs.height,
    additional_nonexact_connectivity_pairs=connectivity_pairs.filter(
        ~pl.col("exact_identity")
    ).height,
    assays=a.height,
    assays_at_least_3=a.filter(pl.col("compounds") >= 3).height,
    assays_at_least_10=a.filter(pl.col("compounds") >= 10).height,
    purpose="Descriptive evidence only. No calibration parameters applied or selected for OOF.",
    matching=(
        "Primary table exact standardized compound identity; "
        "connectivity-only stereo mismatches separate"
    ),
    tox21="Active inhibitor high-quality replicate median; matches curation potency policy",
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
