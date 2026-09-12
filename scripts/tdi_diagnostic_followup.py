"""Paired uncertainty, threshold sensitivity and fold-mean deployment diagnostics."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["oof"])
old = pl.read_parquet(inp / cfg["old"])
oracle = pl.read_csv(inp / cfg["oracle"])
cuts = pl.read_csv(inp / cfg["cuts"])
source_cfg = json.loads((inp / cfg["ensemble_config"]).read_text())
parts = []
for src in cfg["sources"]:
    t = pl.read_parquet(inp / src["test"])
    if "probability" in t.columns:
        if "prediction" in t.columns:
            t = t.drop("prediction")
        t = t.rename({"probability": "prediction"})
    parts.append(
        t.select("compound_id", "endpoint", "fold", "prediction", "method").with_columns(
            (pl.lit(src["prefix"] + "::") + pl.col("method")).alias("method")
        )
    )
blind = pl.concat(parts, how="vertical_relaxed")
rows, sensitivity, paired, deployment, threshold_transfer = [], [], [], [], []
rng = np.random.default_rng(731)
for ep, method in cfg["selected"].items():
    t = oof.filter((pl.col("endpoint") == ep) & (pl.col("method") == method)).sort("compound_id")
    b = old.filter(pl.col("endpoint") == ep).sort("compound_id")
    assert t["compound_id"].equals(b["compound_id"])
    assert np.array_equal(t["y_true"], b["y_true"])
    y, p, f = [t[col].to_numpy() for col in ["y_true", "prediction", "fold"]]
    bpred = b["prediction"].to_numpy()
    weighted_auc = sum(
        roc_auc_score(y[f == k], p[f == k]) * sum(y[f == k]) * sum(1 - y[f == k]) for k in range(5)
    ) / sum(sum(y[f == k]) * sum(1 - y[f == k]) for k in range(5))
    opt = oracle.filter((pl.col("endpoint") == ep) & (pl.col("method") == "mean::" + method))
    opt = opt.sort("MCC", descending=True).row(0, named=True)
    row = dict(
        endpoint=ep,
        method=method,
        pooled_auc=float(roc_auc_score(y, p)),
        within_fold_pair_weighted_auc=float(weighted_auc),
        oracle_same_oof_MCC=opt["MCC"],
        oracle_threshold=opt["threshold"],
    )
    rows.append(row)
    groups = [np.flatnonzero((f == k) & (y == v)) for k in range(5) for v in [0, 1]]
    differences = []
    for _ in range(1000):
        ix = np.concatenate([rng.choice(g, len(g), replace=True) for g in groups])
        differences.append(roc_auc_score(y[ix], p[ix]) - roc_auc_score(y[ix], bpred[ix]))
    paired.append(
        dict(
            endpoint=ep,
            AUC_new_minus_old=float(roc_auc_score(y, p) - roc_auc_score(y, bpred)),
            conditional_lo=float(np.quantile(differences, 0.025)),
            conditional_hi=float(np.quantile(differences, 0.975)),
        )
    )
    for k in range(5):
        m = f != k
        sensitivity.append(
            dict(
                endpoint=ep,
                excluded_fold=k,
                n=int(m.sum()),
                AUC=float(roc_auc_score(y[m], p[m])),
                MCC=float(matthews_corrcoef(y[m], p[m] >= t["threshold"].to_numpy()[m])),
            )
        )
    bs = blind.filter(
        (pl.col("endpoint") == ep) & pl.col("method").is_in(source_cfg["ensembles"][method])
    )
    member_n = len(source_cfg["ensembles"][method])
    bs = bs.group_by("compound_id", "endpoint", "fold").agg(
        pl.col("prediction").mean(), pl.len().alias("members")
    )
    assert (bs["members"] == member_n).all() and bs.height == 3750
    cc = cuts.filter((pl.col("endpoint") == ep) & (pl.col("method") == method))
    bs = bs.join(
        cc.select("endpoint", "fold", "threshold"), on=["endpoint", "fold"], validate="m:1"
    )
    for fold in range(5):
        bb = bs.filter(pl.col("fold") == fold)
        yy = t.filter(pl.col("fold") == fold)
        cutoff = bb["threshold"][0]
        threshold_transfer.append(
            dict(
                endpoint=ep,
                fold=fold,
                cutoff=cutoff,
                oof_n=len(yy),
                oof_positive=int((yy["prediction"] >= cutoff).sum()),
                test_positive=int((bb["prediction"] >= cutoff).sum()),
                oof_mean=float(yy["prediction"].mean()),
                test_mean=float(bb["prediction"].mean()),
            )
        )
    deployed = bs.group_by("compound_id", "endpoint").agg(
        pl.col("prediction").mean(),
        pl.col("prediction").std().alias("fold_std"),
        (pl.col("prediction") >= pl.col("threshold")).sum().alias("votes"),
    )
    deployed = deployed.with_columns(pl.lit(cc["threshold"].median()).alias("threshold"))
    deployed = deployed.with_columns(
        (pl.col("prediction") >= pl.col("threshold")).alias("submitted"),
        (pl.col("votes") >= 3).alias("majority_diagnostic_only"),
    )
    deployment.append(deployed)
for name, data in [
    ("ranking-summary", rows),
    ("paired-auc-change", paired),
    ("leave-one-fold-out", sensitivity),
    ("threshold-transfer", threshold_transfer),
]:
    pl.DataFrame(data).write_csv(out / f"{name}.csv")
d = pl.concat(deployment)
d.write_parquet(out / "test-fold-diagnostics.parquet")
s = d.group_by("endpoint").agg(
    pl.col("submitted").sum(),
    pl.col("majority_diagnostic_only").sum(),
    (pl.col("submitted") != pl.col("majority_diagnostic_only"))
    .sum()
    .alias("decision_disagreements"),
    pl.col("fold_std").mean(),
)
s.write_csv(out / "deployment-summary.csv")
print(pl.DataFrame(rows).write_csv())
print(pl.DataFrame(paired).write_csv())
print(s.write_csv())
