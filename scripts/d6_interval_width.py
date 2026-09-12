"""Nested, asymmetric interval-width-conditioned residual point decisions for D6."""

import importlib.util
import json
import os
from pathlib import Path

import httpx
import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


emb, labels, a, nested, old, global_oof = [pl.read_parquet(path(i)) for i in range(6)]
spec = importlib.util.spec_from_file_location("official", path(6))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
for key in ["source_blog", "source_code"]:
    r = httpx.get(cfg[key], follow_redirects=True)
    r.raise_for_status()
    (out / f"{key}.txt").write_text(r.text)
ep = "CYP2D6_pIC50_direct_inhibition"
keys = ["compound_id", "endpoint", "fold"]
ref = global_oof.filter(pl.col("method") == "caruana_weight_crossfit").drop("method").sort(keys)
assert ref.height == 6384
h = (
    labels.filter(pl.col("endpoint") == ep)
    .join(a, on="compound_id", validate="1:1")
    .join(emb, on="compound_id", validate="1:1")
    .sort("compound_id")
)
y, lo, hi = [h[c].to_numpy() for c in ["value", "conf_low", "conf_high"]]
width = np.column_stack([np.log1p(np.maximum(y - lo, 0)), np.log1p(np.maximum(hi - y, 0))])
assert np.all(lo <= y) and np.all(y <= hi)
ecols = sorted(c for c in emb.columns if c.startswith("emb_"))
rows, membership, diagnostics = [], [], []


def fit(x, y, seed, kind, weight=None):
    p = dict(
        objective="regression_l1",
        n_estimators=cfg["base_trees"],
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.05,
        reg_lambda=1,
        n_jobs=3,
        random_state=seed,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    if kind == "width":
        p.update(n_estimators=cfg["width_trees"], num_leaves=7, min_child_samples=50)
    if kind == "calibration":
        p.update(
            n_estimators=cfg["calibration_trees"],
            num_leaves=5,
            min_child_samples=cfg["calibration_min_leaf"],
            reg_lambda=5,
        )
    m = lgb.LGBMRegressor(**p)
    m.fit(x, y, sample_weight=weight)
    return m


for fold in range(5):
    va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
    tr = ~va
    aux = h.select("compound_id").join(
        nested.filter(pl.col("outer_fold") == fold),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    x = np.column_stack(
        [h.select(ecols).to_numpy(), aux.select([f"aux_{j}" for j in range(4)]).to_numpy()]
    )
    assert np.isfinite(x).all()
    ti = np.flatnonzero(tr)
    vi = np.flatnonzero(va)
    inner_p = np.zeros(len(ti))
    inner_w = np.zeros((len(ti), 2))
    for inner, (itr, iva) in enumerate(
        GroupKFold(cfg["inner_splits"]).split(x[tr], groups=h["group"].to_numpy()[tr])
    ):
        m = fit(x[ti[itr]], y[ti[itr]], 42 + fold, "base")
        inner_p[iva] = m.predict(x[ti[iva]])
        m.booster_.save_model(str(out / f"fold{fold}-inner{inner}-base.txt"))
        for side in range(2):
            wm = fit(x[ti[itr]], width[ti[itr], side], 42 + fold, "width")
            inner_w[iva, side] = wm.predict(x[ti[iva]])
            wm.booster_.save_model(str(out / f"fold{fold}-inner{inner}-width{side}.txt"))
        membership.append(
            h[ti[iva]]
            .select("compound_id")
            .with_columns(
                pl.lit(fold).alias("outer_fold"), pl.lit(inner).alias("inner_validation_fold")
            )
        )
    base = fit(x[tr], y[tr], 42 + fold, "base")
    raw = base.predict(x[va])
    base.booster_.save_model(str(out / f"fold{fold}-base.txt"))
    valw = np.zeros((len(vi), 2))
    for side in range(2):
        wm = fit(x[tr], width[tr, side], 42 + fold, "width")
        valw[:, side] = wm.predict(x[va])
        wm.booster_.save_model(str(out / f"fold{fold}-width{side}.txt"))
    preds = {"raw": raw}
    # Width features are predictions only, from inner-excluded models for calibration training.
    for mode in ["point_residual", "interval_residual", "width_interval_residual"]:
        cx = inner_p[:, None]
        vx = raw[:, None]
        if mode == "width_interval_residual":
            cx = np.column_stack([cx, inner_w, inner_w.sum(1), inner_w[:, 1] - inner_w[:, 0]])
            vx = np.column_stack([vx, valw, valw.sum(1), valw[:, 1] - valw[:, 0]])
        if mode == "point_residual":
            cy = y[tr] - inner_p
            tx = cx
            weights = None
        else:
            cy = np.r_[lo[tr] - inner_p, hi[tr] - inner_p]
            tx = np.vstack([cx, cx])
            weights = np.full(len(cy), 0.5)
        cm = fit(tx, cy, 42 + fold, "calibration", weights)
        delta = cm.predict(vx)
        cm.booster_.save_model(str(out / f"fold{fold}-{mode}.txt"))
        for w in cfg["fixed_correction_weights"]:
            preds[f"{mode}_{int(w * 100)}pct"] = raw + w * delta
    canonical = (
        h[vi]
        .select(
            "compound_id",
            "endpoint",
            "fold",
            "is_neighbor",
            "conf_low",
            "conf_high",
            pl.col("value").alias("y_true"),
        )
        .with_columns(pl.lit("20260911T104703655337Z-a3923255").alias("split_id"))
    )
    for method, pred in preds.items():
        rows.append(
            canonical.with_columns(pl.lit(method).alias("method"), pl.Series("prediction", pred))
        )
    diagnostics.append(
        h[vi]
        .select("compound_id", "fold")
        .with_columns(
            pl.Series("lower_log_width_pred", valw[:, 0]),
            pl.Series("upper_log_width_pred", valw[:, 1]),
            pl.Series("lower_log_width_true", width[vi, 0]),
            pl.Series("upper_log_width_true", width[vi, 1]),
        )
    )
    print(f"fold {fold} complete", flush=True)
r = pl.concat(rows)
r.write_parquet(out / "d6-oof.parquet")
pl.concat(membership).write_parquet(out / "inner-memberships.parquet")
d = pl.concat(diagnostics)
d.write_parquet(out / "width-oof.parquet")
control = r.filter(pl.col("method") == "raw").sort(keys)
prev = old.filter(pl.col("endpoint") == ep).sort(keys)
assert control.height == prev.height == 1457
for k in keys + ["y_true", "conf_low", "conf_high", "prediction"]:
    assert control[k].equals(prev[k]), k
complete = []
metrics = []
for method in sorted(r["method"].unique()):
    p = r.filter(pl.col("method") == method).select(keys + ["prediction"])
    z = (
        ref.join(p, on=keys, how="left", suffix="_new", validate="1:1")
        .with_columns(pl.coalesce("prediction_new", "prediction").alias("prediction"))
        .drop("prediction_new")
        .with_columns(pl.lit(method).alias("method"))
    )
    assert z.height == 6384 and z.select(keys).unique().height == 6384
    complete.append(z)
complete.append(ref.with_columns(pl.lit("global_weight_crossfit").alias("method")))
allp = pl.concat(complete, how="diagonal_relaxed")
allp.write_parquet(out / "canonical-oof.parquet")
for (method, endpoint), f in allp.group_by("method", "endpoint"):
    for subset in ["all", "neighbors"]:
        g = f if subset == "all" else f.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            v = g if fold == -1 else g.filter(pl.col("fold") == fold)
            yy, pp, ll, uu = [
                v[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]
            ]
            metrics.append(
                dict(
                    method=method,
                    endpoint=endpoint,
                    subset=subset,
                    fold=fold,
                    n=v.height,
                    ST_RAE=float(
                        scorer.rae_soft_threshold_absolute_error(
                            yy, pp, y_true_lower=ll, y_true_upper=uu
                        )
                    ),
                    MAE=float(np.abs(yy - pp).mean()),
                    Spearman=float(spearmanr(yy, pp).statistic),
                )
            )
m = pl.DataFrame(metrics)
m.write_csv(out / "metrics.csv")
print(m.filter((pl.col("endpoint") == ep) & (pl.col("fold") == -1)).write_csv(), flush=True)
summary = dict(
    config=cfg,
    raw_bitwise_equal=True,
    canonical_rows_per_method=6384,
    d6_rows_per_method=1457,
    actual_outer_bounds_used_only_for_scoring=True,
    width_oof_spearman=[
        float(spearmanr(d[f"{s}_log_width_true"], d[f"{s}_log_width_pred"]).statistic)
        for s in ["lower", "upper"]
    ],
    no_lb_constants=True,
    interpretation=(
        "All base and calibration fits exclude outer endpoint labels; geometry-independent "
        "frozen600 plus outer-excluded primary features; "
        "global reference reuses weight-crossfit models"
    ),
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
