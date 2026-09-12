"""Recorded affine calibration screen on reused OOF predictions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


all_oof = pl.read_parquet(path(0))
test = pl.read_parquet(path(1))
methods = [
    "identity",
    "offset_mse",
    "offset_mae",
    "offset_interval",
    "affine_mse",
    "affine_mae",
    "affine_interval",
    "interval_strength25",
    "interval_strength50",
    "interval_strength75",
    "neighbor_interval",
    "neighbor_weighted_interval",
]
methods += [f"slope_{a}" for a in [0.5, 0.75, 1.25, 1.5, 2.0]]


def fit(t, kind):
    p, y, lo, hi = [t[c].to_numpy() for c in ["prediction", "y_true", "conf_low", "conf_high"]]
    if kind == "identity":
        return 1.0, 0.0
    if kind == "offset_mse":
        return 1.0, float(np.mean(y - p))
    if kind == "offset_mae":
        return 1.0, float(np.median(y - p))
    if kind == "affine_mse":
        a = np.clip(np.cov(p, y, ddof=0)[0, 1] / np.var(p), 0.05, 3.0)
        return float(a), float(y.mean() - a * p.mean())
    if kind == "neighbor_interval":
        return fit(t.filter(pl.col("is_neighbor")), "affine_interval")
    if kind.startswith("interval_strength"):
        a, b = fit(t, "affine_interval")
        strength = int(kind.removeprefix("interval_strength")) / 100
        return 1 + strength * (a - 1), strength * b
    w = np.ones(len(t))
    if kind == "neighbor_weighted_interval":
        w = np.where(t["is_neighbor"].to_numpy(), 3.0, 1.0)
    if kind == "affine_mae":
        lo = hi = y
    bounds = (0.05, 3.0)
    if kind == "offset_interval":
        bounds = (1.0, 1.0)
    if kind.startswith("slope_"):
        a = float(kind.split("_")[1])
        bounds = (a, a)
    n = len(p)
    design = sparse.csr_matrix(np.column_stack([p, np.ones(n)]))
    mat = sparse.vstack(
        [sparse.hstack([design, -sparse.eye(n)]), sparse.hstack([-design, -sparse.eye(n)])]
    ).tocsr()
    sol = linprog(
        np.r_[0.0, 0.0, w / w.sum()],
        A_ub=mat,
        b_ub=np.r_[hi, -lo],
        bounds=[bounds, (-10.0, 10.0)] + [(0, None)] * n,
        method="highs",
    )
    assert sol.success, sol.message
    return float(sol.x[0]), float(sol.x[1])


def score(h):
    y, p, lo, hi = [h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
    err = np.maximum(lo - p, 0) + np.maximum(p - hi, 0)
    den = (np.maximum(lo - y.mean(), 0) + np.maximum(y.mean() - hi, 0)).mean()
    return dict(
        n=len(h),
        ST_RAE=float(err.mean() / den),
        MAE=float(np.abs(y - p).mean()),
        RMSE=float(np.mean((y - p) ** 2) ** 0.5),
        bias=float((p - y).mean()),
        Spearman=float(spearmanr(y, p).statistic),
    )


parts, params, tests, rows = [], [], [], []
for base in cfg["bases"]:
    data = all_oof.filter(pl.col("method") == base)
    assert len(data) == 6384
    for (ep,), h in data.group_by("endpoint"):
        h = h.sort("compound_id")
        for kind in methods:
            pp = []
            for fold in range(5):
                train = h.filter(pl.col("fold") != fold)
                val = h.filter(pl.col("fold") == fold)
                assert set(train["compound_id"]).isdisjoint(val["compound_id"])
                a, b = fit(train, kind)
                params.append(
                    dict(
                        base=base,
                        endpoint=ep,
                        calibration=kind,
                        fold=fold,
                        slope=a,
                        intercept=b,
                        n_fit=len(train),
                    )
                )
                v = val.with_columns(
                    (pl.col("prediction") * a + b).alias("prediction"),
                    pl.lit(kind).alias("calibration"),
                )
                pp.append(v)
                rows.append(dict(base=base, endpoint=ep, calibration=kind, fold=fold, **score(v)))
            pred = pl.concat(pp)
            parts.append(pred)
            rows.append(dict(base=base, endpoint=ep, calibration=kind, fold=-1, **score(pred)))
            a, b = fit(h, kind)
            params.append(
                dict(
                    base=base,
                    endpoint=ep,
                    calibration=kind,
                    fold=-1,
                    slope=a,
                    intercept=b,
                    n_fit=len(h),
                )
            )
            if base == "community_direct_mean":
                t = test.filter(pl.col("endpoint") == ep)
                assert len(t) == 750
                tests.append(
                    t.with_columns(
                        (pl.col("prediction") * a + b).alias("prediction"),
                        pl.lit(kind).alias("calibration"),
                    )
                )
pl.concat(parts).write_parquet(out / "calibrated-oof.parquet")
pl.concat(tests).write_parquet(out / "calibrated-test.parquet")
pl.DataFrame(params).write_csv(out / "coefficients.csv")
metrics = pl.DataFrame(rows).sort("base", "endpoint", "calibration", "fold")
metrics.write_csv(out / "metrics.csv")
metrics.filter(pl.col("fold") == -1).group_by("base", "calibration").agg(
    pl.col("ST_RAE").mean().alias("macro_ST_RAE")
).sort("base", "macro_ST_RAE").write_csv(out / "macro.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            methods=methods,
            crossfit=(
                "Calibration excludes assessed fold labels; reused base OOF, not fully nested"
            ),
            test="Native predictions mapped using all OOF labels; no test labels",
            selection="All fixed methods reported. Best-of-screen is exploratory.",
        ),
        indent=2,
    )
)
print(metrics.filter((pl.col("fold") == -1) & pl.col("endpoint").str.contains("2D6")).write_csv())
