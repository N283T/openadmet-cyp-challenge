"""Audit inner group isolation and conditional paired D6 deltas."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(i):
    return pl.read_parquet(next(inp.glob(f"{i:03d}-*")))


p, m, a, labels = map(read, range(4))
ep = "CYP2D6_pIC50_direct_inhibition"
source = labels.filter(pl.col("endpoint") == ep).join(a, on="compound_id", validate="1:1")
for fold in range(5):
    train = source.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
    mm = m.filter(pl.col("outer_fold") == fold).join(a, on="compound_id", validate="1:1")
    assert set(mm["compound_id"]) == set(train["compound_id"])
    assert mm["compound_id"].n_unique() == mm.height
    assert (
        mm.group_by("group")
        .agg(pl.col("inner_validation_fold").n_unique())["inner_validation_fold"]
        .max()
        == 1
    )
keys = ["compound_id", "endpoint", "fold"]
r = p.filter(pl.col("method") == "raw").sort(keys)
for method in p["method"].unique():
    h = p.filter(pl.col("method") == method).sort(keys)
    for col in keys + ["y_true", "conf_low", "conf_high"]:
        assert h[col].equals(r[col]), (method, col)
    assert h.height == 6384 and np.isfinite(h["prediction"].to_numpy()).all()
d = p.filter(pl.col("endpoint") == ep).join(
    a.select("compound_id", "group"), on="compound_id", validate="m:1"
)
ref = d.filter(pl.col("method") == "raw").sort("compound_id")
groups = sorted(ref["group"].unique())
gidx = np.array([groups.index(g) for g in ref["group"]])
rng = np.random.default_rng(42)
counts = np.array(
    [
        np.bincount(rng.integers(len(groups), size=len(groups)), minlength=len(groups))
        for _ in range(1000)
    ]
)
lo, hi, y = [ref[c].to_numpy() for c in ["conf_low", "conf_high", "y_true"]]
weights = counts[:, gidx]
mean = (weights @ y) / weights.sum(1)
den = (
    weights
    * (np.maximum(lo[None, :] - mean[:, None], 0) + np.maximum(mean[:, None] - hi[None, :], 0))
).sum(1)
loss = {}
for method in d["method"].unique():
    q = d.filter(pl.col("method") == method).sort("compound_id")["prediction"].to_numpy()
    loss[method] = np.maximum(lo - q, 0) + np.maximum(q - hi, 0)
rows = []
for method in loss:
    for baseline in [
        "raw",
        "global_weight_crossfit",
        "interval_residual_50pct",
        "interval_residual_100pct",
    ]:
        delta = weights @ (loss[method] - loss[baseline]) / den
        rows.append(
            dict(
                method=method,
                baseline=baseline,
                delta_low=float(np.quantile(delta, 0.025)),
                delta_high=float(np.quantile(delta, 0.975)),
                delta_mean=float(delta.mean()),
            )
        )
pl.DataFrame(rows).write_csv(out / "conditional-bootstrap.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            canonical_rows=6384,
            d6_rows=1457,
            inner_group_disjoint=True,
            outer_label_exclusion=True,
            finite_predictions=True,
            bootstrap=(
                "1000 fixed-group conditional draws; "
                "excludes model refit and method selection uncertainty"
            ),
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows).write_csv())
