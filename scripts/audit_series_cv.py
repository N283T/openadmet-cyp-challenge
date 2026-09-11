"""Audit common series assignments and seed-excluded evaluation subsets."""

import importlib.util
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


a = read("assignments")
seeds = pl.read_csv(next(inp.glob("*-hit-seeds.csv")), schema_overrides={"group": pl.String})
assert a["compound_id"].n_unique() == a.height
assert a.group_by("group").agg(pl.col("fold").n_unique())["fold"].max() == 1
seedids = set(seeds["compound_id"])
labels = read("measurements").filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
f = labels.join(a, on="compound_id", validate="m:1").with_columns(
    pl.col("compound_id").is_in(list(seedids)).alias("is_seed")
)
rows = []
for subset in ["all", "hit_series", "neighbors_only"]:
    g = f if subset == "all" else f.filter(pl.col("hit_origins") != "")
    if subset == "neighbors_only":
        g = g.filter(~pl.col("is_seed"))
    rows.append(
        g.group_by("endpoint", "fold")
        .agg(
            pl.len().alias("n"),
            pl.col("value").mean().alias("mean"),
            pl.col("value").median().alias("median"),
            (pl.col("value") < 4).sum().alias("below4"),
        )
        .with_columns(pl.lit(subset).alias("subset"))
    )
pl.concat(rows).sort("subset", "endpoint", "fold").write_csv(out / "subset-counts.csv")
summary = f.group_by("endpoint").agg(
    pl.len().alias("all_n"),
    (pl.col("hit_origins") != "").sum().alias("hit_series_n"),
    ((pl.col("hit_origins") != "") & ~pl.col("is_seed")).sum().alias("neighbors_n"),
)
summary.write_csv(out / "subset-summary.csv")
if list(inp.glob("*-oof.parquet")):
    oof = read("oof")
    spec = importlib.util.spec_from_file_location(
        "official", next(inp.glob("*-custom_scoring_functions.py"))
    )
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    records = []
    for (strategy, endpoint, model), g in oof.group_by("strategy", "endpoint", "model"):
        assert g["compound_id"].n_unique() == g.height
        for subset in ["all", "hit_series", "neighbors_only"]:
            h = g if subset == "all" else g.filter(pl.col("hit_origins") != "")
            if subset == "neighbors_only":
                h = h.filter(~pl.col("compound_id").is_in(list(seedids)))
            y, p, lo, hi = [
                h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]
            ]
            records.append(
                dict(
                    strategy=strategy,
                    endpoint=endpoint,
                    model=model,
                    subset=subset,
                    n=h.height,
                    ST_RAE=float(
                        scorer.rae_soft_threshold_absolute_error(
                            y, p, y_true_lower=lo, y_true_upper=hi
                        )
                    ),
                    MAE=float(np.abs(y - p).mean()),
                    Spearman=float(spearmanr(y, p).statistic),
                )
            )
    m = pl.DataFrame(records)
    m.write_csv(out / "metrics.csv")
    m.group_by("strategy", "model", "subset").agg(
        pl.col("ST_RAE", "MAE", "Spearman").mean()
    ).write_csv(out / "macro.csv")
write_json(
    out / "summary.json",
    dict(
        compounds=a.height,
        groups=a["group"].n_unique(),
        seed_groups=seeds["group"].n_unique(),
        group_integrity=True,
        coverage=summary.to_dicts(),
    ),
)
print(summary)
