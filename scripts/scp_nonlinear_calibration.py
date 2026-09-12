"""Fixed nonlinear maps cross-fitted across existing OOF folds (development CV)."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, eye, hstack, vstack
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression


def interpolation(x, knots):
    idx = np.clip(np.searchsorted(knots, x, side="right") - 1, 0, len(knots) - 2)
    fraction = np.clip((x - knots[idx]) / (knots[idx + 1] - knots[idx]), 0, 1)
    a = np.zeros((len(x), len(knots)))
    a[np.arange(len(x)), idx] = 1 - fraction
    a[np.arange(len(x)), idx + 1] = fraction
    return a


def interval_map(p, lo, hi, n_knots, penalty):
    knots = np.unique(np.quantile(p, np.linspace(0, 1, n_knots)))
    a = csr_matrix(interpolation(p, knots))
    n, k = a.shape
    # Slack is distance to the acceptable interval; regularizer anchors the map to identity.
    z = csr_matrix((n, k))
    zz = csr_matrix((k, n))
    monotone = csr_matrix(np.eye(k)[:-1] - np.eye(k, k, 1)[:-1])
    aub = vstack(
        [
            hstack([a, -eye(n), z]),
            hstack([-a, -eye(n), z]),
            hstack([eye(k), zz, -eye(k)]),
            hstack([-eye(k), zz, -eye(k)]),
            hstack([monotone, csr_matrix((k - 1, n + k))]),
        ]
    ).tocsr()
    b = np.r_[hi, -lo, knots, -knots, np.zeros(k - 1)]
    result = linprog(
        np.r_[np.zeros(k), np.full(n, 1 / n), np.full(k, penalty / k)],
        A_ub=aub,
        b_ub=b,
        bounds=[(None, None)] * k + [(0, None)] * (n + k),
        method="highs",
        options={"threads": 1},
    )
    assert result.success, result.message
    values = result.x[:k]
    assert np.diff(values).min() >= -1e-7
    return knots, values


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    source = pl.read_parquet(next(inp.glob("000-*")))
    spec = importlib.util.spec_from_file_location("scorer", next(inp.glob("001-*")))
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    outputs, parameters, memberships = [], [], []
    keys = ["compound_id", "endpoint", "fold"]
    for base in cfg["base_methods"]:
        data = source.filter(pl.col("method") == base).sort(keys)
        assert data.height == 6384 and not data.select(keys).is_duplicated().any()
        for (endpoint,), target in data.group_by("endpoint", maintain_order=True):
            for fold in range(5):
                train = target.filter(pl.col("fold") != fold)
                test = target.filter(pl.col("fold") == fold)
                assert set(train["compound_id"]).isdisjoint(set(test["compound_id"]))
                p, y, lo, hi = [
                    train[c].to_numpy() for c in ["prediction", "y_true", "conf_low", "conf_high"]
                ]
                pt = test["prediction"].to_numpy()
                memberships.append(
                    target.select("compound_id", "endpoint", "fold").with_columns(
                        pl.lit(base).alias("base_method"),
                        pl.lit(fold).alias("calibration_fold"),
                        (pl.col("fold") != fold).alias("calibration_training"),
                    )
                )
                maps = {"identity": pt}
                iso = IsotonicRegression(out_of_bounds="clip").fit(p, y)
                maps["isotonic_l2"] = iso.predict(pt)
                parameters.append(
                    dict(
                        base=base,
                        endpoint=endpoint,
                        fold=fold,
                        method="isotonic_l2",
                        knots=iso.X_thresholds_.tolist(),
                        values=iso.y_thresholds_.tolist(),
                    )
                )
                for bins in [5, 9]:
                    q = np.linspace(0, 1, bins)
                    knots = np.unique(np.quantile(p, q))
                    # Fixed quantile transport estimates train marginal label distribution only.
                    levels = np.searchsorted(np.sort(p), knots, side="right") / len(p)
                    vals = np.quantile(y, levels)
                    maps[f"quantile{bins}"] = np.interp(pt, knots, vals)
                    parameters.append(
                        dict(
                            base=base,
                            endpoint=endpoint,
                            fold=fold,
                            method=f"quantile{bins}",
                            knots=knots.tolist(),
                            values=vals.tolist(),
                        )
                    )
                    for penalty in [0.0, 0.05, 0.2]:
                        name = f"interval_l1_k{bins}_reg{penalty:g}"
                        knots, vals = interval_map(p, lo, hi, bins, penalty)
                        maps[name] = np.interp(pt, knots, vals)
                        parameters.append(
                            dict(
                                base=base,
                                endpoint=endpoint,
                                fold=fold,
                                method=name,
                                knots=knots.tolist(),
                                values=vals.tolist(),
                            )
                        )
                knots, vals = interval_map(p, y, y, 5, 0.05)
                maps["point_l1_k5_reg0.05"] = np.interp(pt, knots, vals)
                parameters.append(
                    dict(
                        base=base,
                        endpoint=endpoint,
                        fold=fold,
                        method="point_l1_k5_reg0.05",
                        knots=knots.tolist(),
                        values=vals.tolist(),
                    )
                )
                for name in ["isotonic_l2", "quantile5", "quantile9"]:
                    maps[name + "_half"] = 0.5 * (maps[name] + pt)
                for name, pred in maps.items():
                    assert np.isfinite(pred).all()
                    outputs.append(
                        test.with_columns(
                            pl.Series("prediction", pred),
                            pl.lit(base).alias("base_method"),
                            pl.lit(name).alias("calibration"),
                            pl.lit(base + "__" + name).alias("method"),
                        )
                    )
                print(f"{base} {endpoint} fold{fold} finished", flush=True)
    predictions = pl.concat(outputs)
    predictions.write_parquet(out / "calibrated-oof.parquet")
    pl.concat(memberships).write_parquet(out / "calibration-membership.parquet")
    (out / "maps.json").write_text(json.dumps(parameters))

    def metrics(t):
        y, p, lo, hi = [t[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
        return dict(
            n=len(y),
            ST_RAE=float(
                scorer.rae_soft_threshold_absolute_error(y, p, y_true_lower=lo, y_true_upper=hi)
            ),
            MAE=float(np.abs(y - p).mean()),
            interval_MAE=float(np.maximum(np.maximum(lo - p, p - hi), 0).mean()),
            bias=float(np.mean(p - y)),
            Spearman=float(spearmanr(y, p).statistic),
        )

    rows = []
    for (base, method, endpoint), t in predictions.group_by(
        "base_method", "calibration", "endpoint"
    ):
        rows.append(dict(base_method=base, calibration=method, endpoint=endpoint, **metrics(t)))
    scores = pl.DataFrame(rows).sort("base_method", "endpoint", "calibration")
    scores.write_csv(out / "metrics.csv")
    scores.group_by("base_method", "calibration").agg(
        pl.col("ST_RAE").mean().alias("macro_ST_RAE")
    ).sort("base_method", "macro_ST_RAE").write_csv(out / "macro-metrics.csv")
    rows = []
    for (base, method, endpoint, fold), t in predictions.group_by(
        "base_method", "calibration", "endpoint", "fold"
    ):
        rows.append(
            dict(base_method=base, calibration=method, endpoint=endpoint, fold=fold, **metrics(t))
        )
    pl.DataFrame(rows).write_csv(out / "fold-metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                rows=predictions.height,
                maps=len(parameters),
                evaluation="Cross-fold calibration of reused OOF; development comparison, "
                "not fully nested; all fixed methods reported",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
