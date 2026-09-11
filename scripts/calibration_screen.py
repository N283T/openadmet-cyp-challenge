"""Exploratory second-stage calibration on fixed OOF; NOT fully nested base-model CV."""

import importlib.util
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import write_json


def fit_map(kind, p, y, lo, hi, w):
    if kind == "identity":
        return {"a": 1.0, "b": 0.0}
    if kind == "isotonic":
        m = IsotonicRegression(out_of_bounds="clip").fit(p, y, sample_weight=w)
        return {"x": m.X_thresholds_.tolist(), "z": m.y_thresholds_.tolist()}
    if kind == "binned":
        edges = np.unique(np.quantile(p, [0, 0.2, 0.4, 0.6, 0.8, 1]))
        bins = np.searchsorted(edges[1:-1], p)
        x, z, ww = [], [], []
        for b in np.unique(bins):
            ix = bins == b
            x.append(float(np.average(p[ix], weights=w[ix])))
            z.append(float(np.average(y[ix], weights=w[ix])))
            ww.append(float(w[ix].sum()))
        m = IsotonicRegression().fit(x, z, sample_weight=ww)
        return {"x": list(x), "z": m.predict(x).tolist()}
    if kind == "affine_mse":
        m = LinearRegression(positive=True).fit(p[:, None], y, sample_weight=w)
        return {"a": float(m.coef_[0]), "b": float(m.intercept_)}
    # Exact convex linear program for interval absolute error or ordinary MAE.
    if kind == "affine_mae":
        lo = hi = y
    n = len(p)
    x = sparse.csr_matrix(np.column_stack([p, np.ones(n)]))
    mat = sparse.vstack(
        [sparse.hstack([x, -sparse.eye(n)]), sparse.hstack([-x, -sparse.eye(n)])]
    ).tocsr()
    slope_bounds = (1, 1) if kind == "shift_st" else (0.05, 3.0)
    res = linprog(
        np.r_[0.0, 0.0, w / w.sum()],
        A_ub=mat,
        b_ub=np.r_[hi, -lo],
        bounds=[slope_bounds, (-10.0, 10.0)] + [(0, None)] * n,
        method="highs",
    )
    assert res.success, res.message
    return {"a": float(res.x[0]), "b": float(res.x[1])}


def apply_map(m, p):
    return m["a"] * p + m["b"] if "a" in m else np.interp(p, m["x"], m["z"])


def soft(p, lo, hi):
    return np.maximum(lo - p, 0) + np.maximum(p - hi, 0)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(snapshot_input(inp, name + ".parquet"))

    data = read("regression-oof-diagnostics").sort("endpoint", "compound_id")
    eps = sorted(data["endpoint"].unique())
    test = read("ensemble_probabilities").filter(pl.col("endpoint").is_in(eps))
    comps = read("compounds")
    splits = read("assignments").filter(pl.col("strategy") == "umap")
    spec = importlib.util.spec_from_file_location(
        "official", snapshot_input(inp, "custom_scoring_functions.py")
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    assert data.height == 6525 and test.height == 3000
    assert data.select("compound_id", "endpoint").n_unique() == data.height
    assert set(data["compound_id"]).isdisjoint(test["compound_id"])
    structures = dict(comps.select("compound_id", "standardized_smiles").iter_rows())
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    needed = sorted(set(data["compound_id"]) | set(test["compound_id"]))
    fps, keys = {}, {}
    for cid in needed:
        mol = Chem.MolFromSmiles(structures[cid])
        assert mol is not None
        fps[cid] = gen.GetFingerprintAsNumPy(mol)
        keys[cid] = Chem.MolToInchiKey(mol).split("-")[0]
    cluster = dict(splits.select("compound_id", "group").iter_rows())
    methods = [("identity", False)] + [(k, iw) for k in cfg["maps"] for iw in [False, True]]
    scores, params, weights, domain, predictions, testpred, oracle, scenarios = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    bootstrap_data = {}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for ep, ax in zip(eps, axes.flat, strict=True):
        f = data.filter(pl.col("endpoint") == ep)
        t = test.filter(pl.col("endpoint") == ep).sort("compound_id")
        ids, tids = f["compound_id"].to_list(), t["compound_id"].to_list()
        p, y, lo, hi = [f[c].to_numpy() for c in ["prediction", "y_true", "conf_low", "conf_high"]]
        tp = t["prediction"].to_numpy().astype(float)
        fold = f["fold"].to_numpy()
        group = np.array([cluster[c] for c in ids])
        # Domain-label cross-fitting only; no activity labels enter weight estimation.
        allids = ids + tids
        x = np.array([fps[c] for c in allids])
        d = np.r_[np.zeros(len(ids)), np.ones(len(tids))]
        g = np.array([keys[c] for c in allids])
        prob, ratio = np.zeros(len(d)), np.zeros(len(d))
        for tr, va in StratifiedGroupKFold(5, shuffle=True, random_state=cfg["seed"]).split(
            x, d, g
        ):
            assert set(g[tr]).isdisjoint(g[va])
            model = LogisticRegression(
                C=0.1, solver="liblinear", max_iter=1000, random_state=cfg["seed"]
            )
            model.fit(x[tr], d[tr])
            prob[va] = model.predict_proba(x[va])[:, 1]
            ratio[va] = (
                prob[va] / np.maximum(1 - prob[va], 1e-6) * (d[tr] == 0).sum() / (d[tr] == 1).sum()
            )
        w = np.clip(ratio[: len(ids)], 1 / 3, 3)
        w /= w.mean()
        domain.append(
            dict(
                endpoint=ep,
                auc=float(roc_auc_score(d, prob)),
                ess=float(w.sum() ** 2 / (w * w).sum()),
                n=len(w),
                clipped_low=float((ratio[: len(ids)] < 1 / 3).mean()),
                clipped_high=float((ratio[: len(ids)] > 3).mean()),
            )
        )
        weights.append(
            f.select("compound_id", "endpoint", "fold").with_columns(
                pl.Series("weight", w), pl.Series("domain_probability", prob[: len(ids)])
            )
        )

        def metric(q, ix, y=y, lo=lo, hi=hi):
            den = soft(np.full(ix.sum(), y[ix].mean()), lo[ix], hi[ix]).sum()
            st = float(soft(q[ix], lo[ix], hi[ix]).sum() / den)
            official_st = official.rae_soft_threshold_absolute_error(
                y[ix], q[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
            )
            assert np.isclose(st, official_st)
            return dict(
                ST_RAE=st,
                MAE=float(np.abs(q[ix] - y[ix]).mean()),
                Spearman=float(spearmanr(y[ix], q[ix]).statistic),
                bias=float((q[ix] - y[ix]).mean()),
            )

        qmaps = {}
        for kind, iw in methods:
            name = kind + ("_iw" if iw else "")
            q = np.empty_like(p)
            for held in range(5):
                va = fold == held
                # Guard direct calibration fitting against connectivity duplicates across folds.
                forbidden = {keys[c] for c in np.array(ids)[va]}
                tr = (fold != held) & np.array([keys[c] not in forbidden for c in ids])
                m = fit_map(kind, p[tr], y[tr], lo[tr], hi[tr], w[tr] if iw else np.ones(tr.sum()))
                q[va] = apply_map(m, p[va])
                params.append(
                    dict(endpoint=ep, method=name, fold=held, fit_n=int(tr.sum()), model=m)
                )
            m = fit_map(kind, p, y, lo, hi, w if iw else np.ones(len(p)))
            params.append(dict(endpoint=ep, method=name, fold=-1, fit_n=len(p), model=m))
            fullq = apply_map(m, p)
            oracle.append(
                dict(endpoint=ep, method=name, **metric(fullq, np.ones(len(p), dtype=bool)))
            )
            tq = apply_map(m, tp)
            assert np.isfinite(q).all() and np.isfinite(tq).all()
            for held in [-1, 0, 1, 2, 3, 4]:
                ix = np.ones(len(p), dtype=bool) if held == -1 else fold == held
                scores.append(dict(endpoint=ep, method=name, fold=held, **metric(q, ix)))
            predictions.append(
                f.select(
                    "compound_id", "endpoint", "fold", "split_id", "y_true", "conf_low", "conf_high"
                ).with_columns(pl.Series("prediction", q), pl.lit(name).alias("method"))
            )
            testpred.append(
                t.select("compound_id", "endpoint").with_columns(
                    pl.Series("prediction", tq),
                    pl.lit(name).alias("method"),
                    pl.lit(-1).alias("fold"),
                    pl.lit("test:five-fold-mean").alias("split_id"),
                )
            )
            qmaps[name] = q
        bootstrap_data[ep] = (y, lo, hi, p, qmaps, group)
        # Fixed shift/stretch scenarios: OOF sensitivity, never a claimed test score.
        surface = np.empty((len(cfg["shifts"]), len(cfg["scales"])))
        for i, shift in enumerate(cfg["shifts"]):
            for j, scale in enumerate(cfg["scales"]):
                q = p.mean() + scale * (p - p.mean()) + shift
                tq = tp.mean() + scale * (tp - tp.mean()) + shift
                met = metric(q, np.ones(len(p), dtype=bool))
                surface[i, j] = met["ST_RAE"]
                scenarios.append(
                    dict(
                        endpoint=ep,
                        shift=shift,
                        scale=scale,
                        **met,
                        test_mean=float(tq.mean()),
                        test_sd=float(tq.std()),
                        test_below4=float((tq < 4).mean()),
                    )
                )
        im = ax.imshow(
            surface,
            origin="lower",
            aspect="auto",
            extent=[min(cfg["scales"]), max(cfg["scales"]), min(cfg["shifts"]), max(cfg["shifts"])],
        )
        ax.scatter([1], [0], c="red", marker="x")
        ax.set(
            title=ep.split("_")[0] + ": OOF ST-RAE",
            xlabel="Scale around prediction mean",
            ylabel="Shift in pIC50",
        )
        fig.colorbar(im, ax=ax)
    fig.savefig(out / "affine-sensitivity.png", dpi=150)
    plt.close(fig)
    # Shared cluster resampling; conditional on fixed trained base models and calibrators.
    rng = np.random.default_rng(cfg["seed"])
    groups = sorted(set(g for v in bootstrap_data.values() for g in v[-1]))
    samples = rng.multinomial(
        len(groups), np.ones(len(groups)) / len(groups), size=cfg["bootstrap_draws"]
    )
    boot = {}
    for ep, (y, lo, hi, p, qmaps, group) in bootstrap_data.items():
        ix = np.array([groups.index(g) for g in group])
        ww = samples[:, ix]
        means = ww @ y / ww.sum(axis=1)
        den = (
            ww
            * (
                np.maximum(lo[None, :] - means[:, None], 0)
                + np.maximum(means[:, None] - hi[None, :], 0)
            )
        ).sum(axis=1)
        for name, q in qmaps.items():
            boot[ep, name] = ww @ (soft(q, lo, hi) - soft(p, lo, hi)) / den
    bootrows = []
    for name in qmaps:
        for ep in eps + ["macro"]:
            draws = (
                np.mean([boot[e, name] for e in eps], axis=0) if ep == "macro" else boot[ep, name]
            )
            a, b = np.quantile(draws, [0.025, 0.975])
            bootrows.append(dict(endpoint=ep, method=name, delta_low=float(a), delta_high=float(b)))
    sc = pl.DataFrame(scores)
    macro = (
        sc.filter(pl.col("fold") == -1)
        .group_by("method")
        .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
        .sort("ST_RAE")
    )
    for name, table in [
        ("metrics", sc),
        ("macro", macro),
        ("same_sample_fit_diagnostic", pl.DataFrame(oracle)),
        ("domain_diagnostics", pl.DataFrame(domain)),
        ("sensitivity", pl.DataFrame(scenarios)),
        ("bootstrap", pl.DataFrame(bootrows)),
    ]:
        table.write_csv(out / (name + ".csv"))
    pl.concat(predictions).write_parquet(out / "calibration-oof.parquet")
    pl.concat(testpred).write_parquet(out / "calibration-test.parquet")
    pl.concat(weights).write_parquet(out / "importance-weights.parquet")
    write_json(out / "parameters.json", params)
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            macro=macro.to_dicts(),
            limitations=[
                "Second-stage cross-fitting is NOT nested base-model CV: "
                "other-fold base models saw held-fold labels.",
                "No hidden test labels; scenario test scores cannot be computed.",
                "Identity weights/OOF calibration do not identify conditional label shift.",
                "Comparisons and bootstrap are exploratory; "
                "no correction for selection or retraining uncertainty.",
                "Test uses five-member average; OOF one held-out model.",
            ],
        ),
    )
    print(macro)


if __name__ == "__main__":
    main()
