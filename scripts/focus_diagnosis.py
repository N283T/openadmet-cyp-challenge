"""Saved OOF residual diagnostics; exploratory cross-fold interval calibration."""

import json
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, eye, hstack, vstack
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
files = sorted(inp.iterdir())
keys = ["compound_id", "endpoint", "fold"]
b = pl.read_parquet(files[0]).sort(keys)
methods = cfg["baseline_methods"][:]
for idx, name in enumerate(cfg["added_methods"], 1):
    x = pl.read_parquet(files[idx]).sort(keys)
    assert x.height == b.height == 6384
    assert x.select(keys).equals(b.select(keys)), name
    for c in ["y_true", "conf_low", "conf_high"]:
        assert x[c].equals(b[c]), (name, c)
    b = b.with_columns(x["prediction"].alias(name))
    methods.append(name)
assert not b.select(keys).is_duplicated().any()


def score(h, pred):
    y, lo, hi = (h[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"])
    err = np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)
    den = (np.maximum(lo - y.mean(), 0) + np.maximum(y.mean() - hi, 0)).mean()
    return {
        "n": h.height,
        "st_rae": float(err.mean() / den),
        "mae": float(np.abs(pred - y).mean()),
        "bias": float((pred - y).mean()),
        "y_mean": float(y.mean()),
        "pred_mean": float(pred.mean()),
        "y_sd": float(y.std()),
        "pred_sd": float(pred.std()),
        "spearman": float(spearmanr(y, pred).statistic),
    }


def calibrate(h, affine):
    p = h["final"].to_numpy()
    n = len(p)
    # a*p+b-slack <= upper; -a*p-b-slack <= -lower
    design = csr_matrix(np.column_stack([p, np.ones(n)]))
    mat = vstack([hstack([design, -eye(n)]), hstack([-design, -eye(n)])]).tocsr()
    sol = linprog(
        np.r_[0.0, 0.0, np.ones(n) / n],
        A_ub=mat,
        b_ub=np.r_[h["conf_high"].to_numpy(), -h["conf_low"].to_numpy()],
        bounds=[(0.5, 1.5) if affine else (1.0, 1.0), (-0.5, 0.5)] + [(0, None)] * n,
        method="highs",
    )
    assert sol.success, sol.message
    return sol.x[:2]


metrics = []
bins = []
corr = []
coeff = []
calparts = []
for (ep,), h in b.group_by("endpoint"):
    for subset, g in [
        ("all", h),
        ("neighbor", h.filter(pl.col("is_neighbor"))),
        ("nonneighbor", h.filter(~pl.col("is_neighbor"))),
    ]:
        for fold in [-1, 0, 1, 2, 3, 4]:
            z = g if fold == -1 else g.filter(pl.col("fold") == fold)
            if z.height < 3:
                continue
            for m in methods:
                metrics.append(
                    dict(
                        endpoint=ep, subset=subset, fold=fold, method=m, **score(z, z[m].to_numpy())
                    )
                )
        for m1, m2 in combinations(methods, 2):
            y = g["y_true"].to_numpy()
            p = g[m1].to_numpy()
            q = g[m2].to_numpy()
            corr.append(
                dict(
                    endpoint=ep,
                    subset=subset,
                    method1=m1,
                    method2=m2,
                    pred_corr=float(np.corrcoef(p, q)[0, 1]),
                    residual_corr=float(np.corrcoef(p - y, q - y)[0, 1]),
                )
            )
    for m in methods:
        for bin_on in ["y_true", m]:
            vals = h[bin_on].to_numpy()
            edges = np.unique(np.quantile(vals, [0, 0.2, 0.4, 0.6, 0.8, 1]))
            ids = np.searchsorted(edges[1:-1], vals, side="right")
            z = h.with_columns(pl.Series("bin", ids))
            for (bid,), g in z.group_by("bin"):
                bins.append(
                    dict(
                        endpoint=ep,
                        method=m,
                        bin_on="truth" if bin_on == "y_true" else "prediction",
                        bin=int(bid),
                        bin_low=float(edges[bid]),
                        bin_high=float(edges[bid + 1]),
                        **score(g, g[m].to_numpy()),
                    )
                )
    for fold in range(5):
        train = h.filter(pl.col("fold") != fold)
        test = h.filter(pl.col("fold") == fold)
        for mode in ["offset", "affine"]:
            a, intercept = calibrate(train, mode == "affine")
            coeff.append(
                dict(
                    endpoint=ep,
                    fold=fold,
                    mode=mode,
                    slope=float(a),
                    intercept=float(intercept),
                    n_train=train.height,
                )
            )
            calparts.append(
                test.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]).with_columns(
                    pl.lit(mode).alias("method"),
                    (a * test["final"] + intercept).alias("prediction"),
                )
            )
cal = pl.concat(calparts)
for (ep, mode), h in cal.group_by("endpoint", "method"):
    for subset, g in [
        ("all", h),
        ("neighbor", h.filter(pl.col("is_neighbor"))),
        ("nonneighbor", h.filter(~pl.col("is_neighbor"))),
    ]:
        for fold in [-1, 0, 1, 2, 3, 4]:
            z = g if fold == -1 else g.filter(pl.col("fold") == fold)
            if z.height >= 3:
                metrics.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        method="final_" + mode,
                        **score(z, z["prediction"].to_numpy()),
                    )
                )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
pl.DataFrame(bins).write_csv(out / "bins.csv")
pl.DataFrame(corr).write_csv(out / "correlations.csv")
pl.DataFrame(coeff).write_csv(out / "calibration-coefficients.csv")
cal.write_parquet(out / "calibrated-oof.parquet")
b.write_parquet(out / "aligned-oof.parquet")
(out / "limitations.json").write_text(
    json.dumps(
        {
            "calibration": (
                "Exploratory only. Calibration held-out fold labels excluded, but predictions "
                "on other OOF folds originate from models potentially trained on this held-out "
                "fold. This is not fully nested end-to-end calibration validation."
            ),
            "rows": b.height,
            "split_changed": False,
        },
        indent=2,
    )
)
print(
    pl.DataFrame(metrics)
    .filter((pl.col("subset") == "all") & (pl.col("fold") == -1))
    .select("endpoint", "method", "st_rae", "bias", "pred_sd", "y_sd")
)
