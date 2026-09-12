"""Verify optional scale priors are invariant to held-out labels and match audited math."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from external_assay_scaling import GROUP_COLUMNS, shrunk_assay_std, training_variance_pools

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(
        next(p for p in inp.iterdir() if p.name.partition("-")[2] == name + ".parquet")
    )


o = read("observations")
a = o.group_by("compound_id", "assay_id").agg(
    *[pl.col(c).first() for c in GROUP_COLUMNS],
    pl.when((pl.col("lower_bound") == pl.col("upper_bound")).fill_null(False).all())
    .then(pl.col("lower_bound").median())
    .otherwise(pl.col("lower_bound").max())
    .alias("lower_bound"),
    pl.when((pl.col("lower_bound") == pl.col("upper_bound")).fill_null(False).all())
    .then(pl.col("upper_bound").median())
    .otherwise(pl.col("upper_bound").min())
    .alias("upper_bound"),
)
a = a.with_columns(
    pl.when(pl.col("lower_bound") > pl.col("upper_bound")).then(None).otherwise(pl.col(c)).alias(c)
    for c in ["lower_bound", "upper_bound"]
)
training = read("membership").filter(pl.col("role") == "training")["compound_id"].to_list()
base = training_variance_pools(a, training)
is_training = pl.col("compound_id").is_in(pl.Series(training).implode())
perturbed = a.with_columns(
    pl.when(~is_training).then(pl.col(c) * (-1000) + 100000).otherwise(pl.col(c)).alias(c)
    for c in ["lower_bound", "upper_bound"]
)
heldout = training_variance_pools(perturbed, training)
assert base.equals(heldout)
# Verify nonvacuity: perturb training potencies in one assay and see its pool change.
assay = (
    a.filter(is_training & (pl.col("lower_bound") == pl.col("upper_bound")))
    .group_by("assay_id")
    .len()
    .sort("len", descending=True)["assay_id"][0]
)
changed = a.with_columns(
    pl.when(is_training & (pl.col("assay_id") == assay))
    .then(pl.col(c) * 2)
    .otherwise(pl.col(c))
    .alias(c)
    for c in ["lower_bound", "upper_bound"]
)
assert not base.equals(training_variance_pools(changed, training))
assert np.isclose(shrunk_assay_std([4.0, 6.0], 0.5, 20, 0.3), np.sqrt(12 / 22))
assert shrunk_assay_std([5.0, 5.0], 0.0, 20, 0.3) == 0.3
assert shrunk_assay_std([], 1.0, 20, 0.3) == 1.0
# Prior audit independently derived the within-assay grouped variances.
expected = read("assay-scale-audit")
pool_lookup = {
    tuple(r[k] for k in GROUP_COLUMNS): r["pool_variance"] for r in base.iter_rows(named=True)
}
for r in expected.iter_rows(named=True):
    assert np.isclose(pool_lookup[tuple(r[k] for k in GROUP_COLUMNS)], r["pooled_within_std"] ** 2)
base.write_parquet(out / "verified-training-variance-pools.parquet")
result = dict(
    passed=True,
    pools=base.height,
    training_compounds=len(training),
    audited_assays=expected.height,
    validation_perturbation=(
        "All held-out lower/upper labels multiply -1000 then add100000; "
        "priors exactly unchanged"
    ),
    training_perturbation_changes_pool=True,
    formula_case_std=shrunk_assay_std([4.0, 6.0], 0.5, 20, 0.3),
    minimum_sigma=0.3,
    default_policy="assay unchanged",
)
(out / "summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
