"""Summarize train generalization gap and paired-assay restricted size contrast."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, log_loss, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
original = pl.read_parquet(inp / cfg["original"])
paired = pl.read_parquet(inp / cfg["paired"])
data = pl.concat([original, paired])
rows = []
for (method, ep, scope, fold), t in data.group_by("method", "endpoint", "scope", "outer_fold"):
    y, p = t["y_true"].to_numpy(), t["prediction"].to_numpy()
    rows.append(
        dict(
            method=method,
            endpoint=ep,
            scope=scope,
            fold=fold,
            n=len(y),
            positive_rate=float(y.mean()),
            AUROC=roc_auc_score(y, p),
            AP=average_precision_score(y, p),
            log_loss=log_loss(y, p),
        )
    )
detail = pl.DataFrame(rows)
detail.write_csv(out / "fold-generalization.csv")
detail.group_by("method", "endpoint", "scope").agg(
    pl.col("n").mean().alias("mean_n"),
    pl.col("AUROC").mean(),
    pl.col("AP").mean(),
    pl.col("log_loss").mean(),
).write_csv(out / "mean-fold-generalization.csv")
pooled = []
for (method, ep), t in data.filter(pl.col("scope") == "validation").group_by("method", "endpoint"):
    pooled.append(
        dict(
            method=method,
            endpoint=ep,
            n=t.height,
            positives=int(t["y_true"].sum()),
            positive_rate=float(t["y_true"].mean()),
            AUROC=roc_auc_score(t["y_true"], t["prediction"]),
            AP=average_precision_score(t["y_true"], t["prediction"]),
            log_loss=log_loss(t["y_true"], t["prediction"]),
        )
    )
pl.DataFrame(pooled).write_csv(out / "pooled-validation.csv")
selected = pl.read_parquet(inp / cfg["selected"])
route = pl.read_parquet(inp / cfg["route"])
allowed = set(
    route.filter(
        (pl.col("endpoint") == "CYP3A4_is_TDI")
        & (pl.col("availability") == "both")
        & pl.col("y_true").is_finite()
    )["compound_id"]
)
assert len(allowed) == 2334 and set(selected["compound_id"]).issubset(allowed)
assert not selected.filter(pl.col("fold") == pl.col("outer_fold")).height
base = paired.filter(
    (pl.col("scope") == "validation") & (pl.col("method") == "paired_assay_full")
).sort("compound_id")
assert base.height == 2334 and set(base["compound_id"]) == allowed
others = []
for seed in [17, 42, 73]:
    t = paired.filter(
        (pl.col("scope") == "validation") & (pl.col("method") == f"paired_assay_d6_size_seed{seed}")
    ).sort("compound_id")
    assert t.select("compound_id", "y_true", "fold").equals(
        base.select("compound_id", "y_true", "fold")
    )
    others.append(t["prediction"].to_numpy())
    for fold in range(5):
        n = original.filter(
            (pl.col("method") == "full")
            & (pl.col("endpoint") == "CYP2D6_is_TDI")
            & (pl.col("scope") == "train")
            & (pl.col("outer_fold") == fold)
        ).height
        assert (
            selected.filter(
                (pl.col("method") == f"paired_assay_d6_size_seed{seed}")
                & (pl.col("outer_fold") == fold)
            ).height
            == n
        )
y, p = base["y_true"].to_numpy(), base["prediction"].to_numpy()
strata = []
for fold in range(5):
    strata.append(
        [
            np.flatnonzero((base["fold"].to_numpy() == fold) & (base["group"].to_numpy() == g))
            for g in base.filter(pl.col("fold") == fold)["group"].unique()
        ]
    )
rng = np.random.default_rng(2026)
replicates = []
for rep in range(500):
    ix = np.concatenate(
        [
            np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
            for groups in strata
        ]
    )
    auc = roc_auc_score(y[ix], p[ix]) - np.mean([roc_auc_score(y[ix], q[ix]) for q in others])
    ap = average_precision_score(y[ix], p[ix]) - np.mean(
        [average_precision_score(y[ix], q[ix]) for q in others]
    )
    replicates.append(dict(replicate=rep, AUROC_difference=auc, AP_difference=ap))
pl.DataFrame(replicates).write_csv(out / "paired-assay-bootstrap-replicates.csv")
cirows = []
for metric, func in [("AUROC", roc_auc_score), ("AP", average_precision_score)]:
    point = func(y, p) - np.mean([func(y, q) for q in others])
    ci = np.quantile([r[f"{metric}_difference"] for r in replicates], [0.025, 0.975])
    cirows.append(
        dict(
            comparison="A4_both_assays_full_minus_D6size",
            metric=metric,
            difference=point,
            ci_low=ci[0],
            ci_high=ci[1],
            n=2334,
            clusters=sum(map(len, strata)),
        )
    )
pl.DataFrame(cirows).write_csv(out / "paired-assay-bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            paired_assay_membership_exact=True,
            fixed_validation_2334=True,
            train_logloss_included=True,
            no_new_models=True,
        ),
        indent=2,
    )
)
