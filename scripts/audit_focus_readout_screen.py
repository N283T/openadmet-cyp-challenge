"""Independently verify focused OOF and group-bootstrap fixed blend differences."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(name):
    return next(inp.glob("*-" + name))


new = pl.read_parquet(path("focus-oof.parquet"))
base = pl.read_parquet(path("external-oof.parquet"))
a = pl.read_parquet(path("parent-assignments.parquet"))
spec = importlib.util.spec_from_file_location("official", path("custom_scoring_functions.py"))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
keys = ["compound_id", "endpoint", "fold"]
for (ep, method), data in new.group_by("endpoint", "method"):
    h = data.sort(keys)
    ref = base.filter(pl.col("endpoint") == ep).sort(keys)
    assert h.height == ref.height
    for col in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert h[col].equals(ref[col]), col
    h = h.join(
        a.select("compound_id", "group"), on="compound_id", validate="1:1", maintain_order="left"
    )
    y, lo, hi, p0, p1 = (
        h["y_true"].to_numpy(),
        h["conf_low"].to_numpy(),
        h["conf_high"].to_numpy(),
        ref["prediction"].to_numpy(),
        h["prediction"].to_numpy(),
    )
    _, groups = np.unique(h["group"].to_numpy(), return_inverse=True)
    count = np.bincount(groups)
    ysum = np.bincount(groups, weights=y)
    error0 = np.maximum(lo - p0, 0) + np.maximum(p0 - hi, 0)
    for weight in [0.25, 0.5]:
        pred = (1 - weight) * p0 + weight * p1
        actual = scorer.rae_soft_threshold_absolute_error(y, pred, y_true_lower=lo, y_true_upper=hi)
        original = scorer.rae_soft_threshold_absolute_error(y, p0, y_true_lower=lo, y_true_upper=hi)
        delta_error = np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0) - error0
        diff_by_group = np.bincount(groups, weights=delta_error)
        rng = np.random.default_rng(42)
        draws = []
        for _ in range(1000):
            weights = np.bincount(rng.integers(0, len(count), len(count)), minlength=len(count))
            n = np.dot(weights, count)
            mean = np.dot(weights, ysum) / n
            denominator = np.dot(
                weights[groups], np.maximum(lo - mean, 0) + np.maximum(mean - hi, 0)
            )
            draws.append(np.dot(weights, diff_by_group) / denominator)
        rows.append(
            dict(
                endpoint=ep,
                method=method,
                new_weight=weight,
                n=h.height,
                baseline=original,
                score=actual,
                delta=actual - original,
                delta_ci_low=float(np.quantile(draws, 0.025)),
                delta_ci_high=float(np.quantile(draws, 0.975)),
                fraction_improved=float(np.mean(np.array(draws) < 0)),
            )
        )
pl.DataFrame(rows).write_csv(out / "verified-comparison.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            exact_keys_labels_bounds=True,
            official_scorer=True,
            bootstrap_draws=1000,
            qualification=(
                "Conditional group bootstrap of fixed development predictions; "
                "not selection-adjusted"
            ),
            matched_per_method=2829,
            comparisons=rows,
        ),
        indent=2,
    )
)
print(json.dumps(rows, indent=2))
