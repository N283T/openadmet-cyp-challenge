"""Matched point versus credible-interval kernel regression for CYP2D6."""

import importlib.util
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


spec = importlib.util.spec_from_file_location("official", path(5))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
threadpool_limits(3)
kernels = np.load(path(0))
ids, labels, assignments, baseline = [pl.read_parquet(path(i)) for i in range(1, 5)]
ep = "CYP2D6_pIC50_direct_inhibition"
h = ids.join(
    labels.filter(pl.col("endpoint") == ep), on="compound_id", validate="1:1", maintain_order="left"
).join(assignments, on="compound_id", validate="1:1", maintain_order="left")
assert ids["compound_id"].equals(h["compound_id"])
y, lo, hi = [h[c].to_numpy() for c in ["value", "conf_low", "conf_high"]]
assert np.isfinite(np.column_stack([y, lo, hi])).all()
groups, anchors = h["group"].to_numpy(), h["is_anchor"].to_numpy()


def score(frame, p):
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["value"].to_numpy(),
            p,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


def fit(k, lower, upper, reg):
    # Squared distance to the measured interval, with an RKHS norm penalty.
    n = len(lower)

    factor = cholesky(k + np.eye(n) * 1e-9, lower=True)

    def fun(v):
        w, b = v[:-1], v[-1]
        p = factor @ w + b
        d = np.minimum(p - lower, 0) + np.maximum(p - upper, 0)
        loss = np.mean(d * d) / 2 + reg * np.dot(w, w) / 2
        grad = np.r_[factor.T @ d / n + reg * w, d.mean()]
        return loss, grad

    v0 = np.r_[np.zeros(n), np.median((lower + upper) / 2)]
    result = minimize(
        fun,
        v0,
        jac=True,
        method="L-BFGS-B",
        options=dict(maxiter=1000, ftol=1e-11, gtol=1e-6, maxcor=20),
    )
    assert np.isfinite(result.fun) and np.isfinite(result.x).all()
    assert result.success, result.message
    result.x[:-1] = solve_triangular(factor.T, result.x[:-1], lower=False)
    return result


rows, trials, fitted, memberships = [], [], [], []
for fold in range(5):
    mask = (h["fold"].to_numpy() == fold) & ~anchors
    tr, va = np.flatnonzero(~mask), np.flatnonzero(mask)
    eligible = sorted(set(groups[tr]) - set(groups[tr[anchors[tr]]]))
    ig = set(
        np.random.default_rng(42 + fold).permutation(eligible)[: max(1, int(0.2 * len(eligible)))]
    )
    it, iv = (
        np.array([i for i in tr if groups[i] not in ig]),
        np.array([i for i in tr if groups[i] in ig]),
    )
    memberships.append(
        h.select("compound_id").with_columns(
            pl.lit(fold).alias("outer_fold"),
            pl.Series("outer_validation", mask),
            pl.Series("inner_validation", np.isin(np.arange(len(h)), iv)),
        )
    )
    for name in cfg["kernels"]:
        k = kernels[name]
        for loss in ["point", "interval"]:
            lower, upper = (y, y) if loss == "point" else (lo, hi)
            scored = []
            for reg in cfg["regularization"]:
                result = fit(k[np.ix_(it, it)], lower[it], upper[it], reg)
                pred = k[np.ix_(iv, it)] @ result.x[:-1] + result.x[-1]
                val = score(h[iv], pred)
                trials.append(
                    dict(
                        fold=fold,
                        method=name + "_" + loss,
                        reg=reg,
                        ST_RAE=val,
                        converged=bool(result.success),
                        iterations=result.nit,
                    )
                )
                scored.append((val, reg))
            _, reg = min(scored)
            result = fit(k[np.ix_(tr, tr)], lower[tr], upper[tr], reg)
            p = k[np.ix_(va, tr)] @ result.x[:-1] + result.x[-1]
            method = name + "_" + loss
            joblib.dump(
                dict(
                    alpha=result.x[:-1],
                    intercept=result.x[-1],
                    reference_ids=h[tr]["compound_id"].to_list(),
                    reg=reg,
                ),
                out / f"{method}-{fold}.joblib",
            )
            fitted.append(
                dict(
                    method=method,
                    fold=fold,
                    reg=reg,
                    converged=bool(result.success),
                    iterations=result.nit,
                    message=str(result.message),
                )
            )
            rows.append(
                h[va]
                .select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                )
                .with_columns(
                    pl.lit(method).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", p),
                )
            )
            print(fitted[-1], score(h[va], p), flush=True)
oof = pl.concat(rows)
oof.write_parquet(out / "d6-interval-kernel-oof.parquet")
pl.DataFrame(trials).write_csv(out / "inner-trials.csv")
pl.concat(memberships).write_parquet(out / "memberships.parquet")
base = baseline.filter(
    (pl.col("endpoint") == ep) & (pl.col("method") == "caruana_weight_crossfit")
).sort("compound_id")
metrics = []
for method in oof["method"].unique().sort():
    p = oof.filter(pl.col("method") == method).sort("compound_id")
    for col in ["compound_id", "fold", "y_true", "conf_low", "conf_high"]:
        assert p[col].equals(base[col]), col
    for w in [0, 0.25, 0.5, 1]:
        pred = (1 - w) * base["prediction"].to_numpy() + w * p["prediction"].to_numpy()
        metrics.append(
            dict(
                method=method,
                weight=w,
                ST_RAE=score(p.rename({"y_true": "value"}), pred),
                Spearman=float(spearmanr(p["y_true"], pred).statistic),
            )
        )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(json.dumps(dict(config=cfg, fits=fitted), indent=2))
