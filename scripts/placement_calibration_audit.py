"""Independent affine reconstruction and conditional paired group bootstrap."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


raw, affine, nonlinear = [pl.read_parquet(path(i)) for i in range(3)]
coef, metrics = pl.read_csv(path(3)), pl.read_csv(path(4))
labels, assigned = pl.read_parquet(path(5)), pl.read_parquet(path(6))
spec = importlib.util.spec_from_file_location("official", path(7))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
base = raw.filter(pl.col("method") == "community_direct_mean").sort(keys)
expected = (
    labels.filter(pl.col("endpoint").is_in(base["endpoint"].unique()))
    .join(assigned, on="compound_id", validate="m:1")
    .filter(~pl.col("is_anchor"))
    .sort(keys)
)
assert len(expected) == len(base) == 6384
for c in keys:
    assert expected[c].equals(base[c])
for a, b in [("value", "y_true"), ("conf_low", "conf_low"), ("conf_high", "conf_high")]:
    np.testing.assert_allclose(expected[a], base[b], rtol=0, atol=1e-9)
checks = []
for (method, endpoint, calibration), h in affine.group_by("method", "endpoint", "calibration"):
    source = raw.filter((pl.col("method") == method) & (pl.col("endpoint") == endpoint))
    for fold in range(5):
        hh = h.filter(pl.col("fold") == fold).sort(keys)
        pp = source.filter(pl.col("fold") == fold).sort(keys)
        cc = coef.filter(
            (pl.col("base") == method)
            & (pl.col("endpoint") == endpoint)
            & (pl.col("calibration") == calibration)
            & (pl.col("fold") == fold)
        )
        assert cc.height == 1
        for c in keys + ["y_true", "conf_low", "conf_high"]:
            assert hh[c].equals(pp[c])
        assert cc["n_fit"][0] == source.filter(pl.col("fold") != fold).height
        np.testing.assert_allclose(
            hh["prediction"],
            pp["prediction"] * cc["slope"][0] + cc["intercept"][0],
            atol=1e-12,
            rtol=0,
        )
    for fold in [-1, 0, 1, 2, 3, 4]:
        hh = h if fold == -1 else h.filter(pl.col("fold") == fold)
        y, p, lo, hi = [hh[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
        score = float(
            scorer.rae_soft_threshold_absolute_error(y, p, y_true_lower=lo, y_true_upper=hi)
        )
        m = metrics.filter(
            (pl.col("base") == method)
            & (pl.col("endpoint") == endpoint)
            & (pl.col("calibration") == calibration)
            & (pl.col("fold") == fold)
        )
        np.testing.assert_allclose(score, m["ST_RAE"][0], rtol=0, atol=1e-12)
        checks.append(
            dict(base=method, endpoint=endpoint, calibration=calibration, fold=fold, ST_RAE=score)
        )
pl.DataFrame(checks).write_csv(out / "recomputed-metrics.csv")
d6 = (
    base.filter(pl.col("endpoint").str.contains("2D6"))
    .join(assigned.select("compound_id", "group"), on="compound_id", validate="1:1")
    .sort(keys)
)
assert len(d6) == 1457
preds = {"identity": d6["prediction"].to_numpy()}
for name, table, basecol in [
    ("offset_interval", affine, "method"),
    ("interval_l1_k9_reg0.05", nonlinear, "base_method"),
]:
    h = table.filter(
        (pl.col(basecol) == "community_direct_mean")
        & (pl.col("calibration") == name)
        & pl.col("endpoint").str.contains("2D6")
    ).sort(keys)
    for c in keys + ["y_true", "conf_low", "conf_high"]:
        assert h[c].equals(d6[c])
    preds[name] = h["prediction"].to_numpy()
y, lo, hi = [d6[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"]]


def scores(ix):
    return {
        k: float(
            scorer.rae_soft_threshold_absolute_error(
                y[ix], p[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
            )
        )
        for k, p in preds.items()
    }


point = scores(np.arange(len(d6)))
groups = d6["group"].to_numpy()
unique = np.unique(groups)
blocks = [np.flatnonzero(groups == g) for g in unique]
rng = np.random.default_rng(20260912)
rows = []
for replicate in range(500):
    selected = rng.integers(0, len(blocks), len(blocks))
    ix = np.concatenate([blocks[i] for i in selected])
    row = dict(replicate=replicate, n=len(ix), **scores(ix))
    rows.append(row)
boot = pl.DataFrame(rows)
boot.write_parquet(out / "paired-bootstrap.parquet")
summary = []
for method in list(preds)[1:]:
    delta = boot[method].to_numpy() - boot["identity"].to_numpy()
    summary.append(
        dict(
            method=method,
            identity_ST_RAE=point["identity"],
            ST_RAE=point[method],
            delta=point[method] - point["identity"],
            delta_ci_low=float(np.quantile(delta, 0.025)),
            delta_ci_high=float(np.quantile(delta, 0.975)),
            fraction_delta_negative=float(np.mean(delta < 0)),
        )
    )
pl.DataFrame(summary).write_csv(out / "bootstrap-summary.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            canonical_rows=6384,
            d6_rows=1457,
            parent_groups=len(blocks),
            bootstrap_replicates=500,
            affine_metric_checks=len(checks),
            prediction_reconstruction=True,
            scorer_exact=True,
            interpretation="Conditional paired group bootstrap of saved development OOF; "
            "maps fixed; methods screened; "
            "not confirmatory and does not include fitting or method-selection uncertainty",
        ),
        indent=2,
    )
)
print(pl.DataFrame(summary).write_csv())
