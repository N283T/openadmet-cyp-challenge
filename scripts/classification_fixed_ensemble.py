"""Fixed probability means with thresholds selected from isolated inner fits."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
banks = {"oof": {}, "inner": {}, "test": {}}
for source in cfg["sources"]:
    for kind in banks:
        t = pl.read_parquet(inp / source[kind])
        if kind == "inner":
            if "fold" in t.columns:
                t = t.drop("fold")
            t = t.rename({"outer_fold": "fold"})
        if "probability" in t.columns:
            if "prediction" in t.columns:
                t = t.drop("prediction")
            t = t.rename({"probability": "prediction"})
        for (method,), h in t.group_by("method"):
            label = source["prefix"] + "::" + method
            cols = keys + ["prediction"] + (["y_true"] if kind != "test" else [])
            h = h.select(cols).sort(keys)
            h = h.with_columns(pl.col("prediction").cast(pl.Float64), pl.col("fold").cast(pl.Int64))
            if kind != "test":
                h = h.with_columns(pl.col("y_true").cast(pl.Int64))
            assert not h.select(keys).is_duplicated().any()
            assert len(h) == (7500 if kind == "test" else 5081), (kind, label, len(h))
            assert h["prediction"].is_finite().all()
            assert h["prediction"].is_between(0.0, 1.0).all()
            banks[kind][label] = h


def choose(y, p):
    order = np.argsort(-p, kind="stable")
    yy = y[order]
    pp = p[order]
    k = np.r_[np.flatnonzero(pp[:-1] != pp[1:]) + 1, len(p)]
    tp = np.cumsum(yy)[k - 1]
    fp = k - tp
    fn = y.sum() - tp
    tn = len(y) - y.sum() - fp
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    score = np.divide(tp * tn - fp * fn, den, out=np.zeros_like(den, dtype=float), where=den > 0)
    options = [
        (float(s), -abs(float(t) - 0.5), float(t)) for s, t in zip(score, pp[k - 1], strict=True)
    ]
    options.append((float(matthews_corrcoef(y, p >= 0.5)), 0.0, 0.5))
    return max(options)[2]


results, innerrows, blindrows, metricrows, thresholds = [], [], [], [], []
for name, members in cfg["ensembles"].items():
    tables = {}
    for kind in banks:
        parts = [banks[kind][member] for member in members]
        ref = parts[0]
        for h in parts[1:]:
            assert h.select(keys).equals(ref.select(keys))
            if kind != "test":
                assert np.array_equal(h["y_true"].to_numpy(), ref["y_true"].to_numpy())
        tables[kind] = ref.with_columns(
            pl.Series("prediction", np.mean([p["prediction"].to_numpy() for p in parts], axis=0)),
            pl.lit(name).alias("method"),
        )
    oof, inner, blind = tables["oof"], tables["inner"], tables["test"]
    innerrows.append(inner)
    parts = []
    for (ep,), h in oof.group_by("endpoint"):
        cuts = []
        for fold in range(5):
            train = inner.filter((pl.col("endpoint") == ep) & (pl.col("fold") == fold))
            val = h.filter(pl.col("fold") == fold)
            assert set(train["compound_id"]).isdisjoint(val["compound_id"])
            cut = choose(train["y_true"].to_numpy().astype(int), train["prediction"].to_numpy())
            cuts.append(cut)
            thresholds.append(dict(method=name, endpoint=ep, fold=fold, threshold=cut))
            parts.append(val.with_columns(pl.lit(cut).alias("threshold")))
        t = (
            blind.filter(pl.col("endpoint") == ep)
            .group_by("compound_id", "endpoint")
            .agg(pl.col("prediction").mean())
            .with_columns(
                pl.lit(name).alias("method"),
                pl.lit(float(np.median(cuts))).alias("threshold"),
                pl.lit(-1).alias("fold"),
            )
        )
        blindrows.append(t)
    result = pl.concat(parts)
    results.append(result)
    for (ep,), h in result.group_by("endpoint"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            z = h if fold == -1 else h.filter(pl.col("fold") == fold)
            y = z["y_true"].to_numpy()
            p = z["prediction"].to_numpy()
            for policy, cut in [("fixed_05", 0.5), ("inner_mcc", z["threshold"].to_numpy())]:
                metricrows.append(
                    dict(
                        method=name,
                        endpoint=ep,
                        fold=fold,
                        policy=policy,
                        MCC=float(matthews_corrcoef(y, p >= cut)),
                        AUROC=float(roc_auc_score(y, p)),
                        AP=float(average_precision_score(y, p)),
                        n=len(z),
                    )
                )
pl.concat(results).write_parquet(out / "ensemble-oof.parquet")
pl.concat(innerrows).write_parquet(out / "ensemble-inner.parquet")
pl.concat(blindrows).with_columns(
    (pl.col("prediction") >= pl.col("threshold")).alias("binary")
).write_parquet(out / "ensemble-test.parquet")
pl.DataFrame(metricrows).write_csv(out / "metrics.csv")
pl.DataFrame(thresholds).write_csv(out / "thresholds.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            selection="Fixed means; thresholds from matching inner models; screened comparisons",
            deployment="Fivefold mean probabilities, median inner threshold; no full-label refit",
        ),
        indent=2,
    )
)
