"""Retain two prespecified 2D6 additions as OOF-only candidates, without promotion."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
base = pl.read_parquet(next(inp.glob("*-" + cfg.get("baseline_filename", "external-oof.parquet"))))
new = pl.read_parquet(next(inp.glob("*-focus-oof.parquet"))).filter(
    (pl.col("method") == "combined_lgbm") & (pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition")
)
keys = ["compound_id", "endpoint", "fold"]
merged = base.join(
    new.select(*keys, pl.col("prediction").alias("new_prediction")),
    on=keys,
    how="left",
    validate="1:1",
)
assert merged.height == 6384 and merged["new_prediction"].drop_nulls().len() == 1457
spec = importlib.util.spec_from_file_location(
    "official", next(inp.glob("*-custom_scoring_functions.py"))
)
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
metrics = []
for weight in [0.25, 0.5]:
    tag = f"focus_2d6_combined_{int(weight * 100)}pct"
    candidate = (
        merged.with_columns(
            pl.when(pl.col("new_prediction").is_not_null())
            .then((1 - weight) * pl.col("prediction") + weight * pl.col("new_prediction"))
            .otherwise(pl.col("prediction"))
            .alias("prediction"),
            pl.lit(tag).alias("method"),
        )
        .drop("new_prediction")
        .sort(keys)
    )
    unchanged = candidate.filter(pl.col("endpoint") != "CYP2D6_pIC50_direct_inhibition")
    original = base.filter(pl.col("endpoint") != "CYP2D6_pIC50_direct_inhibition").sort(keys)
    assert unchanged.select(keys + ["prediction"]).equals(original.select(keys + ["prediction"]))
    candidate.write_parquet(out / f"{tag}.parquet")
    scores = []
    for (ep,), g in candidate.group_by("endpoint"):
        val = float(
            scorer.rae_soft_threshold_absolute_error(
                g["y_true"].to_numpy(),
                g["prediction"].to_numpy(),
                y_true_lower=g["conf_low"].to_numpy(),
                y_true_upper=g["conf_high"].to_numpy(),
            )
        )
        metrics.append(dict(method=tag, endpoint=ep, ST_RAE=val))
        scores.append(val)
    metrics.append(dict(method=tag, endpoint="macro", ST_RAE=float(np.mean(scores))))
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "manifest.json").write_text(
    json.dumps(
        dict(
            status="OOF-only development candidates; not promoted or submitted",
            changed_endpoint="CYP2D6_pIC50_direct_inhibition",
            weights=[0.25, 0.5],
            other_predictions_bitwise_unchanged=True,
            metrics=metrics,
        ),
        indent=2,
    )
)
print(json.dumps(metrics))
