"""Reconstruct kernel OOFs and verify fold masks and selected inner parameters."""

import importlib.util
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


oof, members, assignments, base, ids, trials = [
    pl.read_parquet(path(i)) if i != 5 else pl.read_csv(path(i)) for i in range(6)
]
kernels = np.load(path(6))
spec = importlib.util.spec_from_file_location("official", path(7))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
order = {v: i for i, v in enumerate(ids["compound_id"])}
rows = []
maxdiff = 0.0
for (method,), h in oof.group_by("method"):
    h = h.sort("compound_id")
    b = base.filter(
        (pl.col("method") == "caruana_weight_crossfit") & (pl.col("endpoint") == h["endpoint"][0])
    ).sort("compound_id")
    assert len(h) == len(b) == 1457
    for c in ["compound_id", "endpoint", "fold", "conf_low", "conf_high", "y_true", "is_neighbor"]:
        assert h[c].equals(b[c]), c
    for fold in range(5):
        m = members.filter(pl.col("outer_fold") == fold).join(
            assignments, on="compound_id", validate="1:1"
        )
        expected = (m["fold"].to_numpy() == fold) & ~m["is_anchor"].to_numpy()
        assert np.array_equal(m["outer_validation"], expected)
        training = m.filter(~pl.col("outer_validation"))
        assert not m.filter(pl.col("inner_validation"))["is_anchor"].any()
        ig = set(m.filter(pl.col("inner_validation"))["group"])
        assert ig.isdisjoint(training.filter(~pl.col("inner_validation"))["group"])
        eligible = sorted(
            set(training["group"]) - set(training.filter(pl.col("is_anchor"))["group"])
        )
        expected_ig = set(
            np.random.default_rng(42 + fold).permutation(eligible)[
                : max(1, int(0.2 * len(eligible)))
            ]
        )
        assert ig == expected_ig
        asset = next(p for p in inp.iterdir() if p.name.endswith(f"-{method}-{fold}.joblib"))
        fit = joblib.load(asset)
        assert set(fit["reference_ids"]) == set(training["compound_id"])
        part = h.filter(pl.col("fold") == fold)
        tr = [order[i] for i in fit["reference_ids"]]
        va = [order[i] for i in part["compound_id"]]
        kk = method if cfg["kind"] == "svr" else method.rsplit("_", 1)[0]
        cross = kernels[kk][np.ix_(va, tr)]
        winner = trials.filter((pl.col("method") == method) & (pl.col("fold") == fold))
        if cfg["kind"] == "svr":
            pred = fit["model"].predict(cross)
            best = winner.sort(["ST_RAE", "C", "epsilon"]).row(0, named=True)
            assert fit["model"].C == best["C"] and fit["model"].epsilon == best["epsilon"]
        else:
            pred = cross @ fit["alpha"] + fit["intercept"]
            best = winner.sort(["ST_RAE", "reg"]).row(0, named=True)
            assert fit["reg"] == best["reg"] and winner["converged"].all()
        maxdiff = max(maxdiff, float(np.max(np.abs(pred - part["prediction"].to_numpy()))))
    for w in [0, 0.25, 0.5, 1]:
        p = (1 - w) * b["prediction"].to_numpy() + w * h["prediction"].to_numpy()
        s = float(
            scorer.rae_soft_threshold_absolute_error(
                h["y_true"].to_numpy(),
                p,
                y_true_lower=h["conf_low"].to_numpy(),
                y_true_upper=h["conf_high"].to_numpy(),
            )
        )
        rows.append(
            dict(
                method=method,
                weight=w,
                ST_RAE=s,
                Spearman=float(spearmanr(h["y_true"], p).statistic),
            )
        )
assert maxdiff < 1e-9, maxdiff
pl.DataFrame(rows).write_csv(out / "verified-metrics.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            max_abs_reconstruction=maxdiff,
            canonical_keys=1457,
            exact_memberships=True,
            exact_inner_selection=True,
            rows=rows,
        ),
        indent=2,
    )
)
