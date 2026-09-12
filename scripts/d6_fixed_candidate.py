"""Preserve the prespecified half-strength scope mean and audit its fixed OOF gain."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


old, mean, a = [pl.read_parquet(path(i)) for i in range(3)]
spec = importlib.util.spec_from_file_location("official", path(3))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
old = old.filter(pl.col("method") == "caruana_weight_crossfit").sort(keys)
ep = "CYP2D6_pIC50_direct_inhibition"
b = old.filter(pl.col("endpoint") == ep)
mean = mean.sort(keys)
for c in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert mean[c].equals(b[c]), c
p0, p1 = b["prediction"].to_numpy(), mean["prediction"].to_numpy()
pred = (p0 + p1) / 2
new = pl.concat(
    [old.filter(pl.col("endpoint") != ep), b.with_columns(pl.Series("prediction", pred))]
)
new = new.with_columns(pl.lit("d6_scope_mean_50pct").alias("method")).sort(keys)
assert len(new) == 6384 and not new.select(keys).is_duplicated().any()
new.write_parquet(out / "d6_scope_mean_50pct-oof.parquet")
metrics = []
for method, t in [("previous", old), ("d6_scope_mean_50pct", new)]:
    for (endpoint,), h in t.group_by("endpoint"):
        value = float(
            scorer.rae_soft_threshold_absolute_error(
                h["y_true"].to_numpy(),
                h["prediction"].to_numpy(),
                y_true_lower=h["conf_low"].to_numpy(),
                y_true_upper=h["conf_high"].to_numpy(),
            )
        )
        metrics.append(dict(method=method, endpoint=endpoint, ST_RAE=value))
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
b = b.join(
    a.select("compound_id", "group"), on="compound_id", validate="1:1", maintain_order="left"
)
y, lo, hi = [b[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"]]
_, g = np.unique(b["group"], return_inverse=True)
error = (
    np.maximum(lo - pred, 0)
    + np.maximum(pred - hi, 0)
    - np.maximum(lo - p0, 0)
    - np.maximum(p0 - hi, 0)
)
cnt, ys, es = np.bincount(g), np.bincount(g, weights=y), np.bincount(g, weights=error)
rng = np.random.default_rng(42)
draws = []
for _ in range(2000):
    w = np.bincount(rng.integers(0, len(cnt), len(cnt)), minlength=len(cnt))
    center = np.dot(w, ys) / np.dot(w, cnt)
    denom = np.dot(w[g], np.maximum(lo - center, 0) + np.maximum(center - hi, 0))
    draws.append(np.dot(w, es) / denom)
(out / "verification.json").write_text(
    json.dumps(
        dict(
            keys=6384,
            d6_keys=1457,
            recipe=(
                "50% previous global weight-crossfit +50% prespecified "
                "four-scope LGBM mean for2D6 only"
            ),
            conditional_group_bootstrap_95=np.quantile(draws, [0.025, 0.975]).tolist(),
            qualification=(
                "Conditional on fitted development OOFs; "
                "not selection-adjusted or unseen validation"
            ),
            metrics=metrics,
        ),
        indent=2,
    )
)
