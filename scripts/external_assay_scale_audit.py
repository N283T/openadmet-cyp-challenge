"""Audit train-only assay normalization leverage and candidate pooled/shrunk scales."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(name):
    return next(p for p in inp.iterdir() if p.name.partition("-")[2] == name)


def read(name):
    return pl.read_parquet(path(name + ".parquet"))


tasks = pl.DataFrame(json.loads(path("tasks.json").read_text())).filter(
    pl.col("objective") == "continuous"
)
o = read("observations")
# Match the trainer's exact-replicate medians and mixed-interval intersections.
a = o.group_by("compound_id", "assay_id").agg(
    pl.col("source").first(),
    pl.col("isoform").first(),
    pl.col("measurement_type").first(),
    pl.when((pl.col("lower_bound") == pl.col("upper_bound")).fill_null(False).all())
    .then(pl.col("lower_bound").median())
    .otherwise(pl.col("lower_bound").max())
    .alias("lower"),
    pl.when((pl.col("lower_bound") == pl.col("upper_bound")).fill_null(False).all())
    .then(pl.col("upper_bound").median())
    .otherwise(pl.col("upper_bound").min())
    .alias("upper"),
)
a = a.filter(~(pl.col("lower") > pl.col("upper")).fill_null(False))
a = a.join(read("membership"), on="compound_id", validate="m:1")
a = a.with_columns(
    ((pl.col("lower") == pl.col("upper")) & pl.col("lower").is_finite())
    .fill_null(False)
    .alias("exact")
)
train = a.filter((pl.col("role") == "training") & pl.col("exact"))
groups = ["source", "isoform", "measurement_type"]
pools = train.group_by(groups).agg(
    pl.col("lower").std(ddof=0).alias("pooled_std"),
    pl.len().alias("pooled_n"),
    pl.col("lower").min().alias("pooled_min"),
    pl.col("lower").max().alias("pooled_max"),
)
# Residual variance excludes assay mean shifts from the scale prior.
train = train.with_columns(
    (pl.col("lower") - pl.col("lower").mean().over("assay_id")).alias("centered")
)
within = train.group_by(groups).agg(
    (pl.col("centered") ** 2).mean().sqrt().alias("pooled_within_std")
)
pools = pools.join(within, on=groups, validate="1:1")
pools.write_csv(out / "train-only-pooled-scales.csv")
rows = []
for task in tasks.iter_rows(named=True):
    g = a.filter(pl.col("assay_id") == task["assay_id"])
    tr = g.filter((pl.col("role") == "training") & pl.col("exact"))
    va = g.filter(
        (pl.col("role") == "validation")
        & (pl.col("lower").is_not_null() | pl.col("upper").is_not_null())
    )
    x = tr["lower"].to_numpy()
    rawstd = float(np.std(x))
    assert np.isclose(np.mean(x), task["mean"])
    assert np.isclose(max(0.1, rawstd), task["std"])
    group = g.select(groups).unique().row(0, named=True)
    pool = pools.filter(pl.all_horizontal([pl.col(k) == v for k, v in group.items()])).row(
        0, named=True
    )
    poolstd = max(0.3, pool["pooled_within_std"])
    # Pseudocount20 prior; train-only metadata, no validation values used in scale.
    shrunk = max(0.3, float(np.sqrt((len(x) * rawstd**2 + 20 * poolstd**2) / (len(x) + 20))))
    lo = va["lower"].fill_null(-float("inf")).to_numpy()
    hi = va["upper"].fill_null(float("inf")).to_numpy()
    err = np.maximum(lo - task["mean"], 0) ** 2 + np.maximum(task["mean"] - hi, 0) ** 2
    rows.append(
        task
        | group
        | dict(
            train_exact_n=len(x),
            raw_std=rawstd,
            train_min=float(x.min()),
            train_max=float(x.max()),
            train_range=float(np.ptp(x)),
            inverse_variance=1 / task["std"] ** 2,
            pooled_std=pool["pooled_std"],
            pooled_within_std=poolstd,
            shrunk_std=shrunk,
            shrunk_inverse_variance=1 / shrunk**2,
            validation_known_n=len(err),
            zerohead_raw_interval_mse=float(err.mean()) if len(err) else None,
            zerohead_current_scaled_mse=float(err.mean() / task["std"] ** 2) if len(err) else None,
            zerohead_pooled_scaled_mse=float(err.mean() / poolstd**2) if len(err) else None,
            zerohead_shrunk_scaled_mse=float(err.mean() / shrunk**2) if len(err) else None,
        )
    )
t = pl.DataFrame(rows, infer_schema_length=None)
t.sort("inverse_variance", descending=True).write_parquet(out / "assay-scale-audit.parquet")
t.sort("zerohead_current_scaled_mse", descending=True, nulls_last=True).head(40).write_csv(
    out / "high-leverage-assays.csv"
)
summary = []
for (source,), g in t.group_by("source"):
    summary.append(
        dict(
            source=source,
            assays=g.height,
            train_n_median=g["train_n"].median(),
            train_exact_n_median=g["train_exact_n"].median(),
            train_exact_n_le5=(g["train_exact_n"] <= 5).sum(),
            std_floor_count=(g["std"] <= 0.1000001).sum(),
            std_le025=(g["std"] <= 0.25).sum(),
            std_median=g["std"].median(),
            std_q10=g["std"].quantile(0.1),
            std_q90=g["std"].quantile(0.9),
            inverse_variance_median=g["inverse_variance"].median(),
            inverse_variance_q90=g["inverse_variance"].quantile(0.9),
            zerohead_current_scaled_mse=g["zerohead_current_scaled_mse"].mean(),
            zerohead_raw_interval_mse=g["zerohead_raw_interval_mse"].mean(),
            zerohead_pooled_scaled_mse=g["zerohead_pooled_scaled_mse"].mean(),
            zerohead_shrunk_scaled_mse=g["zerohead_shrunk_scaled_mse"].mean(),
            validation_assays=g["zerohead_current_scaled_mse"].is_not_null().sum(),
        )
    )
pl.DataFrame(summary).write_csv(out / "source-scale-summary.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            source_stats=summary,
            assays=t.height,
            policy=(
                "Descriptive scale leverage audit. Train-only pools and pseudocount20 variance "
                "shrinkage examples; no model trained or predictions changed."
            ),
            scale_groups=groups,
            prior="Within-assay residual pooled variance; minimum sigma0.3; pseudocount20",
        ),
        indent=2,
    )
)
print(json.dumps(summary, indent=2))
