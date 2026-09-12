"""Fit fixed exploratory calibration maps on all OOF and apply to native test means."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
paths = [next(inp.glob(f"{i:03d}-*")) for i in range(4)]
oof = pl.read_parquet(paths[0]).filter(pl.col("method") == "community_direct_mean")
test = pl.read_parquet(paths[1]).with_columns(pl.col("prediction").cast(pl.Float64))
spec = importlib.util.spec_from_file_location("nonlinear", paths[2])
calibrator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(calibrator)
canonical = pl.read_parquet(paths[3])
assert oof.height == 6384 and test.height == 3000
assert test["compound_id"].n_unique() == 750
assert set(test["compound_id"]) == set(canonical["compound_id"])
outputs, maps = [], []
for (endpoint,), t in test.group_by("endpoint"):
    train = oof.filter(pl.col("endpoint") == endpoint)
    assert t.height == 750 and train.height > 1000
    p, lo, hi = [train[c].to_numpy() for c in ["prediction", "conf_low", "conf_high"]]
    pt = t["prediction"].to_numpy()
    outputs.append(t.with_columns(pl.lit("identity").alias("calibration")))
    for name, k, penalty in [("interval_l1_k5_reg0", 5, 0), ("interval_l1_k9_reg0.05", 9, 0.05)]:
        knots, values = calibrator.interval_map(p, lo, hi, k, penalty)
        predicted = np.interp(pt, knots, values)
        assert np.isfinite(predicted).all()
        outputs.append(
            t.with_columns(pl.Series("prediction", predicted), pl.lit(name).alias("calibration"))
        )
        maps.append(
            dict(
                endpoint=endpoint,
                calibration=name,
                knots=knots.tolist(),
                values=values.tolist(),
                train_rows=train.height,
            )
        )
result = pl.concat(outputs).with_columns(
    pl.lit("community_direct_mean").alias("base_method"),
    pl.concat_str(pl.lit("community_direct_mean__"), pl.col("calibration")).alias("method"),
)
assert (
    result.height == 9000
    and not result.select("compound_id", "endpoint", "method").is_duplicated().any()
)
result.write_parquet(out / "calibrated-test.parquet")
(out / "maps.json").write_text(json.dumps(maps, indent=2))
result.group_by("calibration", "endpoint").agg(
    pl.len().alias("n"),
    pl.col("prediction").mean(),
    pl.col("prediction").std().alias("std"),
    pl.col("prediction").min().alias("min"),
    pl.col("prediction").max().alias("max"),
).write_csv(out / "distribution.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            rows=result.height,
            rows_per_method=3000,
            test_compounds=750,
            maps=len(maps),
            methods=["identity", "interval_l1_k5_reg0", "interval_l1_k9_reg0.05"],
            training="All native OOF predictions and available confidence bounds; "
            "test labels unused",
            selection="Fixed exploratory candidates requested after development comparison, "
            "not independently selected",
            submitted=False,
        ),
        indent=2,
    )
)
