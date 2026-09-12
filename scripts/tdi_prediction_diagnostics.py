"""Independent probability diagnostics; no fitting or deployment selection."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
canonical = (
    pl.read_parquet(inp / cfg["canonical"])
    .select("compound_id", "endpoint", "fold", "y_true")
    .with_columns(pl.col("y_true").cast(pl.Int64), pl.col("fold").cast(pl.Int64))
)
keys = ["compound_id", "endpoint", "fold"]
banks = {k: {} for k in ["oof", "inner", "test"]}
for source in cfg["sources"]:
    for kind in banks:
        if kind not in source:
            continue
        t = pl.read_parquet(inp / source[kind])
        t = t.filter(pl.col("endpoint").is_in(canonical["endpoint"].unique()))
        if kind == "inner":
            if "fold" in t.columns:
                t = t.drop("fold")
            t = t.rename({"outer_fold": "fold"})
        if "probability" in t.columns:
            if "prediction" in t.columns:
                t = t.drop("prediction")
            t = t.rename({"probability": "prediction"})
        if "method" not in t.columns:
            t = t.with_columns(pl.lit(source.get("method", "historical")).alias("method"))
        for (method,), h in t.group_by("method"):
            h = h.with_columns(
                pl.col("prediction").cast(pl.Float64), pl.col("fold").cast(pl.Int64)
            ).sort(keys)
            if kind != "test":
                h = h.with_columns(pl.col("y_true").cast(pl.Int64))
            assert not h.select(keys).is_duplicated().any()
            assert h["prediction"].is_between(0, 1).all()
            if kind == "oof":
                assert h.select(canonical.columns).equals(canonical.sort(keys))
            banks[kind][source["prefix"] + "::" + method] = h


def stats(y, p):
    q = np.clip(p, 1e-8, 1 - 1e-8)
    return dict(
        n=len(y),
        positives=int(y.sum()),
        prevalence=float(y.mean()),
        AUROC=float(roc_auc_score(y, p)),
        AP=float(average_precision_score(y, p)),
        Brier=float(np.mean((y - p) ** 2)),
        logloss=float(-np.mean(y * np.log(q) + (1 - y) * np.log1p(-q))),
        pred_mean=float(p.mean()),
        pred_std=float(p.std()),
        MCC05=float(matthews_corrcoef(y, p >= 0.5)),
    )


def distribution(p):
    return {
        "mean": float(p.mean()),
        "std": float(p.std()),
        **{f"q{q}": float(np.quantile(p, q / 100)) for q in [0, 5, 25, 50, 75, 95, 100]},
    }


metrics, thresholds, distributions, reliability, cis, errors, correlations = (
    [],
    [],
    [],
    [],
    [],
    [],
    [],
)
rng = np.random.default_rng(cfg["seed"])
for method, h in banks["oof"].items():
    for (ep,), t in h.group_by("endpoint"):
        y, p, f = t["y_true"].to_numpy(), t["prediction"].to_numpy(), t["fold"].to_numpy()
        base = dict(method=method, endpoint=ep)
        for fold in [-1, 0, 1, 2, 3, 4]:
            mask = np.ones(len(t), bool) if fold == -1 else f == fold
            yy, pp = y[mask], p[mask]
            row = dict(**base, fold=fold, **stats(yy, pp))
            if "threshold" in t.columns:
                cut = t["threshold"].to_numpy()[mask]
                row["MCC_saved_cut"] = float(matthews_corrcoef(yy, pp >= cut))
                row["pred_positive_rate_saved_cut"] = float((pp >= cut).mean())
            # Honest intercept-only comparison: each validation fold uses otherfoldprevalence.
            prior = np.array([y[f != ff].mean() for ff in f[mask]])
            row["Brier_trainprior"] = float(np.mean((yy - prior) ** 2))
            row["logloss_trainprior"] = float(
                -np.mean(yy * np.log(prior) + (1 - yy) * np.log1p(-prior))
            )
            metrics.append(row)
            for cls in [-1, 0, 1]:
                arr = pp if cls == -1 else pp[yy == cls]
                distributions.append(
                    dict(**base, fold=fold, kind="oof", label=cls, **distribution(arr))
                )
        order = np.argsort(p, kind="stable")
        for decile, ix in enumerate(np.array_split(order, 10)):
            reliability.append(
                dict(
                    **base,
                    decile=decile,
                    n=len(ix),
                    pred_mean=float(p[ix].mean()),
                    positive_rate=float(y[ix].mean()),
                )
            )
        if method in cfg["focus"]:
            # Descriptive all-OOF threshold sweep is an optimistic bound, never deployed.
            for cut in np.unique(np.r_[p, 0.5]):
                thresholds.append(
                    dict(
                        **base,
                        kind="oracle_same_oof_only",
                        threshold=float(cut),
                        MCC=float(matthews_corrcoef(y, p >= cut)),
                    )
                )
            samplegroups = [
                np.flatnonzero((f == fold) & (y == cls)) for fold in range(5) for cls in [0, 1]
            ]
            draws = []
            for _ in range(cfg["bootstrap"]):
                ix = np.concatenate([rng.choice(g, len(g), replace=True) for g in samplegroups])
                draws.append(
                    [
                        roc_auc_score(y[ix], p[ix]),
                        matthews_corrcoef(
                            y[ix],
                            p[ix]
                            >= (t["threshold"].to_numpy()[ix] if "threshold" in t.columns else 0.5),
                        ),
                    ]
                )
            draws = np.asarray(draws)
            # Conditional null permutes labels within fold, preserving fold prevalence.
            ranks = rankdata(p)
            null = []
            for _ in range(cfg["permutations"]):
                yp = y.copy()
                for fold in range(5):
                    ix = np.flatnonzero(f == fold)
                    yp[ix] = rng.permutation(yp[ix])
                n1 = yp.sum()
                null.append((ranks[yp == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(y) - n1)))
            cis.append(
                dict(
                    **base,
                    AUROC=float(roc_auc_score(y, p)),
                    AUC_lo=float(np.quantile(draws[:, 0], 0.025)),
                    AUC_hi=float(np.quantile(draws[:, 0], 0.975)),
                    MCC_lo=float(np.quantile(draws[:, 1], 0.025)),
                    MCC_hi=float(np.quantile(draws[:, 1], 0.975)),
                    withinfold_null_mean=float(np.mean(null)),
                    permutation_p=float(
                        (1 + np.sum(np.array(null) >= roc_auc_score(y, p))) / (len(null) + 1)
                    ),
                )
            )
            errors.append(
                t.with_columns(
                    pl.lit(method).alias("method"),
                    (pl.col("prediction") - pl.col("y_true")).abs().alias("abs_error"),
                )
            )
        if method in banks["inner"]:
            inner = banks["inner"][method].filter(pl.col("endpoint") == ep)
            for (fold,), z in inner.group_by("fold"):
                distributions.append(
                    dict(
                        **base,
                        fold=fold,
                        kind="inner",
                        label=-1,
                        **distribution(z["prediction"].to_numpy()),
                    )
                )
        if method in banks["test"]:
            blind = banks["test"][method].filter(pl.col("endpoint") == ep)
            for (fold,), z in blind.group_by("fold"):
                distributions.append(
                    dict(
                        **base,
                        fold=fold,
                        kind="test",
                        label=-1,
                        **distribution(z["prediction"].to_numpy()),
                    )
                )
            avg = blind.group_by("compound_id").agg(pl.col("prediction").mean())
            distributions.append(
                dict(
                    **base,
                    fold=-1,
                    kind="test_mean",
                    label=-1,
                    **distribution(avg["prediction"].to_numpy()),
                )
            )

for ep in canonical["endpoint"].unique():
    names = list(banks["oof"])
    for i, name in enumerate(names):
        a = banks["oof"][name].filter(pl.col("endpoint") == ep).sort(keys)
        for other in names[i + 1 :]:
            b = banks["oof"][other].filter(pl.col("endpoint") == ep).sort(keys)
            assert a.select(keys).equals(b.select(keys))
            correlations.append(
                dict(
                    endpoint=ep,
                    left=name,
                    right=other,
                    pearson=float(np.corrcoef(a["prediction"], b["prediction"])[0, 1]),
                    spearman=float(spearmanr(a["prediction"], b["prediction"]).statistic),
                )
            )
for name, rows in [
    ("metrics", metrics),
    ("distributions", distributions),
    ("reliability", reliability),
    ("conditional_intervals", cis),
    ("oracle_thresholds", thresholds),
    ("correlations", correlations),
]:
    pl.DataFrame(rows).write_csv(out / f"{name}.csv")
pl.concat(
    [
        e.select("compound_id", "endpoint", "fold", "method", "y_true", "prediction", "abs_error")
        for e in errors
    ]
).write_parquet(out / "focus-errors.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            methods=len(banks["oof"]),
            config=cfg,
            uncertainty="Conditional saved fits; stratified row bootstrap, no refits",
            oracle="All-OOF threshold results are fitting-only diagnostics and were not deployed",
        ),
        indent=2,
    )
)
print(pl.DataFrame(cis).write_csv(), flush=True)
