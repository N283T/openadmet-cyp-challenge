"""Audit which training label routes carry AP/control decision swaps."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.classification_evaluation import (
    fixed_rate_calls,
    prepare_predictions,
)

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
assignments = pl.read_parquet(inputs / cfg["assignments"])
route = (
    pl.read_parquet(inputs / cfg["assay_route"])
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", "availability", "shift_regime", "bounds_regime")
)
if route.select("compound_id").is_duplicated().any():
    raise ValueError("Duplicate assay-route compound")

all_swaps = []
for strategy, files in cfg["predictions"].items():
    bank = pl.concat([pl.read_parquet(inputs / name) for name in files]).filter(
        pl.col("split_id") == strategy
    )
    models = cfg["models"][strategy]
    aligned = prepare_predictions(bank, assignments, models, strategy)
    parts = [aligned.filter(pl.col("model") == model).sort("compound_id") for model in models]
    if parts[0].height != cfg["expected_pool_rows"]:
        raise ValueError("Incomplete full pool")
    if parts[0]["compound_id"].to_list() != parts[1]["compound_id"].to_list():
        raise ValueError("Model IDs differ")
    folds = parts[0]["fold"].to_numpy()
    calls = [fixed_rate_calls(part["score"].to_numpy(), folds, cfg["rate"]) for part in parts]
    swaps = (
        parts[0]
        .select("compound_id", "fold", "y_true")
        .with_columns(
            pl.Series("control_call", calls[0]),
            pl.Series("ap_call", calls[1]),
            pl.lit(strategy).alias("strategy"),
        )
    )
    swaps = swaps.with_columns(
        pl.when(pl.col("control_call") & pl.col("ap_call"))
        .then(pl.lit("both"))
        .when(pl.col("control_call"))
        .then(pl.lit("control_only"))
        .when(pl.col("ap_call"))
        .then(pl.lit("ap_only"))
        .otherwise(pl.lit("neither"))
        .alias("decision")
    ).join(route, on="compound_id", how="left", validate="1:1")
    if swaps["availability"].null_count():
        raise ValueError("Assay-route coverage incomplete")
    all_swaps.append(swaps)

result = pl.concat(all_swaps)
result.write_parquet(outputs / "decisions-by-compound.parquet")
summary = (
    result.group_by("strategy", "decision", "availability", "shift_regime", "bounds_regime")
    .agg(
        pl.len().alias("pool_rows"),
        pl.col("y_true").is_not_null().sum().alias("scored_rows"),
        (pl.col("y_true") == 1).sum().alias("positive_labels"),
    )
    .sort("strategy", "decision", "availability", "shift_regime", "bounds_regime")
)
summary.write_csv(outputs / "label-route-summary.csv")
print(summary, flush=True)
