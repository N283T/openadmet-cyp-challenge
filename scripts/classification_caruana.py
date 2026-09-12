"""Nested held-out selection of sparse and Caruana TDI probability ensembles."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
banks = {k: {} for k in ["oof", "inner", "test"]}
for src in cfg["sources"]:
    for kind in banks:
        t = pl.read_parquet(inp / src[kind])
        if kind == "inner":
            if "fold" in t.columns:
                t = t.drop("fold")
            t = t.rename({"outer_fold": "fold"})
        if "probability" in t.columns:
            if "prediction" in t.columns:
                t = t.drop("prediction")
            t = t.rename({"probability": "prediction"})
        for (name,), h in t.group_by("method"):
            h = h.select(*keys, "prediction", *(["y_true"] if kind != "test" else []))
            h = h.with_columns(pl.col("prediction").cast(pl.Float64), pl.col("fold").cast(pl.Int64))
            if kind != "test":
                h = h.with_columns(pl.col("y_true").cast(pl.Int64))
            h = h.sort(keys)
            assert len(h) == (7500 if kind == "test" else 5081)
            assert not h.select(keys).is_duplicated().any()
            assert h["prediction"].is_finite().all() and h["prediction"].is_between(0, 1).all()
            banks[kind][src["prefix"] + "::" + name] = h
names = sorted(banks["oof"])
for kind in banks:
    assert sorted(banks[kind]) == names
    ref = banks[kind][names[0]]
    for h in banks[kind].values():
        assert h.select(keys).equals(ref.select(keys))
        if kind != "test":
            assert h["y_true"].equals(ref["y_true"])


def threshold(y, p):
    order = np.argsort(-p, kind="stable")
    pp, yy = p[order], y[order]
    k = np.r_[np.flatnonzero(pp[:-1] != pp[1:]) + 1, len(p)]
    tp = np.cumsum(yy)[k - 1]
    fp, fn = k - tp, y.sum() - tp
    tn = len(y) - y.sum() - fp
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    scores = np.divide(tp * tn - fp * fn, den, out=np.zeros_like(den, dtype=float), where=den > 0)
    opts = [
        (float(s), -abs(float(t) - 0.5), float(t)) for s, t in zip(scores, pp[k - 1], strict=True)
    ]
    opts.append((float(matthews_corrcoef(y, p >= 0.5)), 0.0, 0.5))
    best = max(opts)
    return best[0], best[2]


def score(y, p, objective):
    if objective == "mcc":
        return threshold(y, p)[0]
    if objective == "brier":
        return -float(np.mean((y - p) ** 2))
    ranks = rankdata(p)
    n = y.sum()
    return float((ranks[y == 1].sum() - n * (n + 1) / 2) / (n * (len(y) - n)))


def fit(x, y, recipe):
    n = x.shape[1]
    best, weights, path = -np.inf, None, []
    if recipe == "fixed_historical":
        weights = np.zeros(n)
        weights[names.index("history::historical_stack_mordred")] = 1.0
        value, cut = threshold(y, x @ weights)
        return weights, cut, dict(objective=value, fit_MCC=value, path=[])
    if recipe == "bagged_sparse3":
        rng = np.random.default_rng(cfg["bagging_seed"])
        groups = [np.flatnonzero(y == cls) for cls in [0, 1]]
        samples = []
        for sample in range(cfg["bagging_samples"]):
            ix = np.concatenate([rng.choice(g, len(g), replace=True) for g in groups])
            w, _, detail = fit(x[ix], y[ix], "caruana_mcc_sparse3")
            samples.append(w)
            path.append(dict(sample=sample, weights=w.tolist(), fit_MCC=detail["fit_MCC"]))
        mean = np.mean(samples, axis=0)
        keep = np.argsort(-mean, kind="stable")[:3]
        weights = np.zeros(n)
        weights[keep] = mean[keep]
        weights /= weights.sum()
        value, cut = threshold(y, x @ weights)
        return weights, cut, dict(objective=value, fit_MCC=value, path=path)
    if recipe.startswith("subset"):
        count = int(recipe[-1])
        for members in itertools.combinations(range(n), count):
            p = x[:, members].mean(axis=1)
            value = score(y, p, "mcc")
            if value > best + 1e-12:
                best = value
                weights = np.zeros(n)
                weights[list(members)] = 1 / count
        path = [dict(step=count, objective=best, members=np.flatnonzero(weights).tolist())]
    else:
        objective = recipe.split("_")[1]
        cap = 3 if recipe.endswith("sparse3") else n
        counts, total = np.zeros(n, int), np.zeros(len(y))
        for step in range(1, cfg["steps"] + 1):
            candidates = [i for i in range(n) if counts[i] or np.count_nonzero(counts) < cap]
            values = [score(y, (total + x[:, i]) / step, objective) for i in candidates]
            choice = candidates[int(np.argmax(values))]
            total += x[:, choice]
            counts[choice] += 1
            value = max(values)
            path.append(dict(step=step, objective=value, added=names[choice]))
            # Earliest strictly best prefix: no automatic preference for 50-member solutions.
            if value > best + 1e-12:
                best, weights = value, counts.copy() / step
    p = x @ weights
    mcc, cut = threshold(y, p)
    return weights, cut, dict(objective=best, fit_MCC=mcc, path=path)


def metric(y, p, binary):
    return dict(
        MCC=float(matthews_corrcoef(y, binary)),
        AUROC=float(roc_auc_score(y, p)),
        AP=float(average_precision_score(y, p)),
        Brier=float(np.mean((y - p) ** 2)),
        positives=int(binary.sum()),
        n=len(y),
    )


rows, tests, deployed, selections, metrics, fitting = [], [], [], [], [], []
for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
    refs, matrices = {}, {}
    for kind in banks:
        tables = [banks[kind][name].filter(pl.col("endpoint") == ep) for name in names]
        refs[kind] = tables[0]
        matrices[kind] = np.column_stack([h["prediction"].to_numpy() for h in tables])
    for recipe in cfg["recipes"]:
        for fold in range(5):
            masks = {kind: refs[kind]["fold"].to_numpy() == fold for kind in banks}
            inner, outer = (
                refs["inner"].filter(pl.Series(masks["inner"])),
                refs["oof"].filter(pl.Series(masks["oof"])),
            )
            assert set(inner["compound_id"]).isdisjoint(outer["compound_id"])
            w, cut, detail = fit(
                matrices["inner"][masks["inner"]], inner["y_true"].to_numpy(), recipe
            )
            selections.append(
                dict(
                    endpoint=ep,
                    recipe=recipe,
                    fold=fold,
                    threshold=cut,
                    weights=dict(zip(names, w.tolist(), strict=True)),
                    **detail,
                )
            )
            p = matrices["oof"][masks["oof"]] @ w
            rows.append(
                outer.select(keys + ["y_true"]).with_columns(
                    pl.lit(recipe).alias("method"),
                    pl.Series("prediction", p),
                    pl.lit(cut).alias("threshold"),
                )
            )
            tests.append(
                refs["test"]
                .filter(pl.Series(masks["test"]))
                .select(keys)
                .with_columns(
                    pl.lit(recipe).alias("method"),
                    pl.Series("prediction", matrices["test"][masks["test"]] @ w),
                    pl.lit(cut).alias("threshold"),
                )
            )
        # Final candidate parameters use all OOF labels. Their fitting score is NOT CV.
        w, cut, detail = fit(matrices["oof"], refs["oof"]["y_true"].to_numpy(), recipe)
        selections.append(
            dict(
                endpoint=ep,
                recipe=recipe,
                fold=-1,
                threshold=cut,
                weights=dict(zip(names, w.tolist(), strict=True)),
                **detail,
            )
        )
        fitting.append(
            dict(
                endpoint=ep,
                method=recipe,
                **metric(
                    refs["oof"]["y_true"].to_numpy(),
                    matrices["oof"] @ w,
                    matrices["oof"] @ w >= cut,
                ),
            )
        )
        blind = (
            refs["test"].select(keys).with_columns(pl.Series("prediction", matrices["test"] @ w))
        )
        blind = (
            blind.group_by("compound_id", "endpoint")
            .agg(pl.col("prediction").mean())
            .with_columns(
                pl.lit(recipe).alias("method"),
                pl.lit(-1).alias("fold"),
                pl.lit(cut).alias("threshold"),
            )
        )
        deployed.append(blind.with_columns((pl.col("prediction") >= cut).alias("binary")))
        print(ep, recipe, "finished", flush=True)
result = pl.concat(rows)
for (ep, method), h in result.group_by("endpoint", "method"):
    for fold in [-1, 0, 1, 2, 3, 4]:
        z = h if fold == -1 else h.filter(pl.col("fold") == fold)
        metrics.append(
            dict(
                endpoint=ep,
                method=method,
                fold=fold,
                **metric(
                    z["y_true"].to_numpy(),
                    z["prediction"].to_numpy(),
                    (z["prediction"] >= z["threshold"]).to_numpy(),
                ),
            )
        )
result.write_parquet(out / "oof.parquet")
pl.concat(tests).write_parquet(out / "outer-weighted-fold-test.parquet")
pl.concat(deployed).write_parquet(out / "deployment.parquet")
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
pl.DataFrame(fitting).write_csv(out / "full-oof-fitting-only.csv")
(out / "selections.json").write_text(json.dumps(selections, indent=2))
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            candidate_names=names,
            evaluation="Weights/cuts from isolated inner fits; outer labels evaluation only",
            deployment="Full OOF-fitted weights/cut on fivefold test means; fitting score separate",
        ),
        indent=2,
    )
)
print(pl.DataFrame(metrics).filter(pl.col("fold") == -1).write_csv())
