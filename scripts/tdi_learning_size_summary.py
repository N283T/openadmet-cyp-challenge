"""Paired cluster bootstrap and membership audit for fixed-validation size study."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
data = pl.read_parquet(inp / cfg["predictions"])
selected = pl.read_parquet(inp / cfg["selected"])
assert not selected.filter(pl.col("fold") == pl.col("outer_fold")).height
assert (
    selected.unique(["compound_id", "endpoint", "method", "outer_fold"]).height == selected.height
)
val = data.filter(pl.col("scope") == "validation")
a4 = val.filter((pl.col("endpoint") == "CYP3A4_is_TDI") & (pl.col("method") == "full")).sort(
    "compound_id"
)
assert a4.height == 3584
for seed in cfg["seeds"]:
    for fold in range(5):
        chosen = selected.filter(
            (pl.col("method") == f"d6_size_seed{seed}") & (pl.col("outer_fold") == fold)
        )
        half = selected.filter(
            (pl.col("method") == f"half_d6_size_seed{seed}") & (pl.col("outer_fold") == fold)
        )
        d6 = selected.filter(
            (pl.col("method") == "full")
            & (pl.col("endpoint") == "CYP2D6_is_TDI")
            & (pl.col("outer_fold") == fold)
        )
        assert chosen.height == d6.height and half.height == d6.height // 2
        assert set(half["compound_id"]).issubset(set(chosen["compound_id"]))
        for method in [f"d6_size_seed{seed}", f"half_d6_size_seed{seed}"]:
            t = val.filter(pl.col("method") == method).sort("compound_id")
            assert t.select("compound_id", "fold", "y_true").equals(
                a4.select("compound_id", "fold", "y_true")
            )
shared = val.filter(pl.col("method") == "shared_population")
for fold in range(5):
    for table, column in [
        (shared, "outer_fold"),
        (selected.filter(pl.col("method") == "shared_population"), "outer_fold"),
    ]:
        p = table.filter((pl.col("endpoint") == "CYP3A4_is_TDI") & (pl.col(column) == fold))
        q = table.filter((pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col(column) == fold))
        assert set(p["compound_id"]) == set(q["compound_id"])


def get(method, ep="CYP3A4_is_TDI"):
    return val.filter((pl.col("method") == method) & (pl.col("endpoint") == ep)).sort("compound_id")


def score(t, idx):
    y, p = t["y_true"].to_numpy()[idx], t["prediction"].to_numpy()[idx]
    return np.array([roc_auc_score(y, p), average_precision_score(y, p)])


comparisons = [
    ("A4_full_minus_D6size", [get("full")] * 3, [get(f"d6_size_seed{s}") for s in cfg["seeds"]]),
    (
        "A4_full_minus_halfD6size",
        [get("full")] * 3,
        [get(f"half_d6_size_seed{s}") for s in cfg["seeds"]],
    ),
    (
        "A4_D6size_minus_halfD6size",
        [get(f"d6_size_seed{s}") for s in cfg["seeds"]],
        [get(f"half_d6_size_seed{s}") for s in cfg["seeds"]],
    ),
    ("shared_A4_minus_D6", [get("shared_population")], [get("shared_population", "CYP2D6_is_TDI")]),
]
rows, bootrows = [], []
for name, left, right in comparisons:
    base = left[0]
    for t in left + right:
        assert t["compound_id"].equals(base["compound_id"])
    strata = []
    for fold in range(5):
        groups = []
        for group in base.filter(pl.col("fold") == fold)["group"].unique():
            groups.append(
                np.flatnonzero(
                    (base["fold"].to_numpy() == fold) & (base["group"].to_numpy() == group)
                )
            )
        strata.append(groups)
    rng = np.random.default_rng(2026)
    samples = []
    for rep in range(cfg["bootstrap_reps"]):
        idx = np.concatenate(
            [
                np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
                for groups in strata
            ]
        )
        if any(len(np.unique(t["y_true"].to_numpy()[idx])) < 2 for t in left + right):
            continue
        delta = np.mean(
            [score(a, idx) - score(b, idx) for a, b in zip(left, right, strict=True)], axis=0
        )
        samples.append(delta)
        bootrows.append(
            dict(comparison=name, replicate=rep, AUROC_difference=delta[0], AP_difference=delta[1])
        )
    point = np.mean(
        [
            score(a, np.arange(base.height)) - score(b, np.arange(base.height))
            for a, b in zip(left, right, strict=True)
        ],
        axis=0,
    )
    for j, metric in enumerate(["AUROC", "AP"]):
        ci = np.quantile(np.array(samples)[:, j], [0.025, 0.975])
        rows.append(
            dict(
                comparison=name,
                metric=metric,
                difference=point[j],
                ci_low=ci[0],
                ci_high=ci[1],
                bootstrap_n=len(samples),
                clusters=sum(map(len, strata)),
                n=base.height,
            )
        )
pl.DataFrame(rows).write_csv(out / "paired-bootstrap.csv")
pl.DataFrame(bootrows).write_parquet(out / "bootstrap-replicates.parquet")
folds = pl.read_csv(inp / cfg["fold_metrics"])
folds.group_by("method", "endpoint", "scope").agg(
    pl.col("n").mean().alias("mean_n"),
    pl.col("AUROC").mean().alias("mean_AUROC"),
    pl.col("AP").mean().alias("mean_AP"),
).write_csv(out / "mean-fold-metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            membership_audit="passed",
            bootstrap="Paired UMAP-group resampling within outer folds",
            seed_summary="Mean metric differences across seeds; no prediction averaging",
            reps=cfg["bootstrap_reps"],
        ),
        indent=2,
    )
)
