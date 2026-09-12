"""Validate saved calibration maps, canonical keys and heldout memberships."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
paths = [next(inp.glob(f"{i:03d}-*")) for i in range(4)]
source, pred, members = map(pl.read_parquet, paths[:3])
parameters = json.loads(paths[3].read_text())
keys = ["compound_id", "endpoint", "fold"]
assert pred.height == 14 * 2 * 6384
errors = []
for (base, calibration), t in pred.group_by("base_method", "calibration"):
    s = source.filter(pl.col("method") == base).sort(keys)
    t = t.sort(keys)
    assert t.select(keys + ["y_true", "conf_low", "conf_high"]).equals(
        s.select(keys + ["y_true", "conf_low", "conf_high"])
    )
    assert t.height == 6384 and np.isfinite(t["prediction"].to_numpy()).all()
    if calibration == "identity":
        assert np.array_equal(t["prediction"].to_numpy(), s["prediction"].to_numpy())
for item in parameters:
    b, e, f, m = [item[c] for c in ["base", "endpoint", "fold", "method"]]
    knots, values = map(np.asarray, [item["knots"], item["values"]])
    assert (np.diff(knots) > 0).all() and np.diff(values).min() >= -1e-7
    s = source.filter(
        (pl.col("method") == b) & (pl.col("endpoint") == e) & (pl.col("fold") == f)
    ).sort(keys)
    t = pred.filter(
        (pl.col("base_method") == b)
        & (pl.col("endpoint") == e)
        & (pl.col("fold") == f)
        & (pl.col("calibration") == m)
    ).sort(keys)
    reconstructed = np.interp(s["prediction"].to_numpy(), knots, values)
    error = float(np.max(np.abs(reconstructed - t["prediction"].to_numpy())))
    assert error < 1e-10
    errors.append(error)
for (_, _, fold), t in members.group_by("base_method", "endpoint", "calibration_fold"):
    assert (t["calibration_training"] == (t["fold"] != fold)).all()
    train = set(t.filter(pl.col("calibration_training"))["compound_id"])
    val = set(t.filter(~pl.col("calibration_training"))["compound_id"])
    assert train and val and train.isdisjoint(val)
(out / "audit.json").write_text(
    json.dumps(
        dict(
            prediction_rows=pred.height,
            fixed_methods=14,
            bases=2,
            maps_checked=len(errors),
            max_reconstruction_error=max(errors),
            memberships="All40 endpoint/fold/base partitions exclude the evaluated fold",
            keys_bounds_identity="Exact canonical preservation",
            status="passed",
        ),
        indent=2,
    )
)
print((out / "audit.json").read_text())
