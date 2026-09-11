"""Check blends and conditional corrections within common base-model holdouts."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from calibration_screen import apply_map, fit_map
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy import sparse
from scipy.optimize import linprog, minimize
from scipy.special import expit
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import write_json


def context(p, chem, kind, state=None):
    if kind == "prediction":
        raw = np.column_stack([p, p * p, np.maximum(p - 4.5, 0), np.maximum(p - 5.5, 0)])
    else:
        raw = np.column_stack([p, chem])
    if state is None:
        state = (raw.mean(axis=0), np.maximum(raw.std(axis=0), 0.01))
    x = np.clip((raw - state[0]) / state[1], -5, 5)
    return np.column_stack([np.ones(len(p)), x]), state


def residual_fit(p, y, lo, hi, x, penalty, mse=False):
    if mse:
        model = Ridge(alpha=penalty * len(y), fit_intercept=False).fit(x, y - p)
        return model.coef_

    # Smooth interval distance for optimization; exact interval loss used in scoring.
    def objective(beta):
        q = p + x @ beta
        z1, z2 = (lo - q) / 0.03, (q - hi) / 0.03
        loss = 0.03 * (np.logaddexp(0, z1) + np.logaddexp(0, z2)).mean()
        return loss + penalty * (beta @ beta), x.T @ (-expit(z1) + expit(z2)) / len(
            p
        ) + 2 * penalty * beta

    result = minimize(
        objective, np.zeros(x.shape[1]), jac=True, method="L-BFGS-B", bounds=[(-1, 1)] * x.shape[1]
    )
    assert result.success, result.message
    return result.x


def convex_fit(x, lo, hi, penalty):
    n, k = x.shape
    # Anchor penalty discourages leaving member zero without interval-loss evidence.
    design = sparse.csr_matrix(np.column_stack([x, np.ones(n)]))
    mat = sparse.vstack(
        [sparse.hstack([design, -sparse.eye(n)]), sparse.hstack([-design, -sparse.eye(n)])]
    ).tocsr()
    obj = np.r_[0, np.full(k - 1, penalty), 0, np.ones(n) / n]
    eq = np.zeros((1, k + 1 + n))
    eq[0, :k] = 1
    result = linprog(
        obj,
        A_ub=mat,
        b_ub=np.r_[hi, -lo],
        A_eq=eq,
        b_eq=[1.0],
        bounds=[(0, 1)] * k + [(-0.5, 0.5)] + [(0, None)] * n,
        method="highs",
    )
    assert result.success, result.message
    return result.x[:k], result.x[k]


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    f = pl.read_parquet(snapshot_input(inp, "regression-oof-diagnostics.parquet")).sort(
        "endpoint", "compound_id"
    )
    eps = sorted(f["endpoint"].unique())
    names = []
    for entry in cfg["members"]:
        name = entry["name"]
        names.append(name)
        pred = pl.read_parquet(inp / entry["input"])
        pred = pred.filter(pl.col("endpoint").is_in(eps))
        assert pred.select("compound_id", "endpoint").n_unique() == 6525
        pred = pred.select(
            "compound_id",
            "endpoint",
            pl.col("prediction").alias(name),
            pl.col("fold").alias("check_fold"),
            pl.col("y_true").alias("check_y"),
        )
        f = f.join(pred, on=["compound_id", "endpoint"], validate="1:1")
        assert f.height == 6525 and (f["fold"] == f["check_fold"]).all()
        assert np.allclose(f["y_true"], f["check_y"])
        f = f.drop("check_fold", "check_y")
    compounds = pl.read_parquet(snapshot_input(inp, "compounds.parquet"))
    structural, keys = {}, {}
    funcs = [
        Descriptors.MolWt,
        Descriptors.MolLogP,
        Descriptors.TPSA,
        Descriptors.NumHDonors,
        Descriptors.NumHAcceptors,
        Descriptors.NumRotatableBonds,
        Descriptors.RingCount,
        Descriptors.FractionCSP3,
    ]
    for cid, smi in compounds.select("compound_id", "standardized_smiles").iter_rows():
        mol = Chem.MolFromSmiles(smi)
        structural[cid] = [fun(mol) for fun in funcs]
        keys[cid] = Chem.MolToInchiKey(mol).split("-")[0]
    assignments = pl.read_parquet(snapshot_input(inp, "assignments.parquet")).filter(
        pl.col("strategy") == "umap"
    )
    groupmap = dict(assignments.select("compound_id", "group").iter_rows())
    spec = importlib.util.spec_from_file_location(
        "official", snapshot_input(inp, "custom_scoring_functions.py")
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    variants = ["identity", "affine_st_half", "equal_top3", "equal_all"] + [
        "blend50_" + name for name in names[1:]
    ]
    variants += [
        f"{kind}_{lam}"
        for kind in ["convex", "prediction_st", "chemical_st", "chemical_mse"]
        for lam in cfg["penalties"]
    ]
    rows, outputs, params, correlations = [], [], [], []
    for ep in eps:
        frame = f.filter(pl.col("endpoint") == ep)
        p, y, lo, hi = [
            frame[c].to_numpy().astype(float)
            for c in ["prediction", "y_true", "conf_low", "conf_high"]
        ]
        mat = frame.select(names).to_numpy().astype(float)
        assert np.allclose(mat[:, 0], p)
        ids = frame["compound_id"].to_numpy()
        folds = frame["fold"].to_numpy()
        groups = np.array([groupmap[c] for c in ids])
        conn = np.array([keys[c] for c in ids])
        chem = np.array([structural[c] for c in ids])
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                correlations.append(
                    dict(
                        endpoint=ep,
                        first=names[i],
                        second=names[j],
                        residual_corr=float(np.corrcoef(mat[:, i] - y, mat[:, j] - y)[0, 1]),
                    )
                )
        for repeat in range(cfg["repeats"]):
            parts = []
            for fold in range(5):
                available = np.unique(groups[folds == fold])
                shuffled = np.random.default_rng(cfg["seed"] + repeat * 100 + fold).permutation(
                    available
                )
                half = set(shuffled[: len(shuffled) // 2])
                side = np.array([g in half for g in groups])
                for s in [False, True]:
                    va = (folds == fold) & (side == s)
                    tr = (folds == fold) & (side != s) & ~np.isin(conn, conn[va])
                    assert tr.sum() > 20 and va.sum() > 20 and set(conn[tr]).isdisjoint(conn[va])
                    parts.append((fold, int(s), tr, va))
            for name in variants:
                if name == "identity":
                    q = p.copy()
                elif name == "equal_top3":
                    q = mat[:, :3].mean(axis=1)
                elif name == "equal_all":
                    q = mat.mean(axis=1)
                elif name.startswith("blend50_"):
                    q = (p + mat[:, names.index(name.removeprefix("blend50_"))]) / 2
                else:
                    q = np.full(len(p), np.nan)
                    for fold, side, tr, va in parts:
                        if name == "affine_st_half":
                            m = fit_map(
                                "affine_st", p[tr], y[tr], lo[tr], hi[tr], np.ones(tr.sum())
                            )
                            q[va] = (p[va] + apply_map(m, p[va])) / 2
                        elif name.startswith("convex"):
                            a, b = convex_fit(mat[tr], lo[tr], hi[tr], float(name.split("_")[-1]))
                            q[va] = mat[va] @ a + b
                            m = dict(weights=a.tolist(), shift=float(b))
                        else:
                            kind = "prediction" if name.startswith("prediction") else "chemical"
                            x, state = context(p[tr], chem[tr], kind)
                            xt, _ = context(p[va], chem[va], kind, state)
                            beta = residual_fit(
                                p[tr],
                                y[tr],
                                lo[tr],
                                hi[tr],
                                x,
                                float(name.split("_")[-1]),
                                mse="mse" in name,
                            )
                            q[va] = p[va] + xt @ beta
                            m = dict(
                                beta=beta.tolist(),
                                center=state[0].tolist(),
                                scale=state[1].tolist(),
                            )
                        params.append(
                            dict(
                                endpoint=ep,
                                method=name,
                                repeat=repeat,
                                fold=fold,
                                half=side,
                                model=m,
                            )
                        )
                assert np.isfinite(q).all()
                for fold in [-1, 0, 1, 2, 3, 4]:
                    ix = np.ones(len(p), dtype=bool) if fold == -1 else folds == fold
                    st = float(
                        official.rae_soft_threshold_absolute_error(
                            y[ix], q[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
                        )
                    )
                    rows.append(
                        dict(
                            endpoint=ep,
                            method=name,
                            repeat=repeat,
                            fold=fold,
                            ST_RAE=st,
                            MAE=float(np.abs(q[ix] - y[ix]).mean()),
                            Spearman=float(spearmanr(y[ix], q[ix]).statistic),
                        )
                    )
                outputs.append(
                    frame.select(
                        "compound_id", "endpoint", "fold", "split_id", "y_true"
                    ).with_columns(
                        pl.Series("prediction", q),
                        pl.lit(name).alias("method"),
                        pl.lit(repeat).alias("repeat"),
                    )
                )
    scores = pl.DataFrame(rows)
    macro = (
        scores.filter(pl.col("fold") == -1)
        .group_by("method", "repeat")
        .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    )
    summary = (
        macro.group_by("method")
        .agg(
            pl.col("ST_RAE", "MAE", "Spearman").mean(),
            pl.col("ST_RAE").min().alias("min"),
            pl.col("ST_RAE").max().alias("max"),
        )
        .sort("ST_RAE")
    )
    scores.write_csv(out / "context-metrics.csv")
    summary.write_csv(out / "context-summary.csv")
    pl.DataFrame(correlations).write_csv(out / "residual-correlations.csv")
    pl.concat(outputs).write_parquet(out / "context-oof.parquet")
    write_json(out / "parameters.json", params)
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            results=summary.to_dicts(),
            limitations=[
                "Within-holdout cluster-half check: small calibration sets; "
                "same-fold representativeness differs from deployment.",
                "Repeated splits reuse labels; method choice is exploratory.",
                "Fixed blends use ordinary UMAP OOF; "
                "fitted corrections use cluster-half calibration.",
            ],
        ),
    )
    print(summary)


if __name__ == "__main__":
    main()
