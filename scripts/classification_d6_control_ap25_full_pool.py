"""Compose the fixed 75/25 control/AP rank blend on both outer splits."""

import json
import os
from pathlib import Path

import polars as pl

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
all_parts = []
for strategy, files in cfg["predictions"].items():
    models = [cfg["control_model"], cfg["candidate_models"][strategy]]
    sources = [pl.read_parquet(inputs / file) for file in files]
    parts = [
        source.filter((pl.col("split_id") == strategy) & (pl.col("model") == model)).sort(
            "compound_id"
        )
        for source, model in zip(sources, models, strict=True)
    ]
    if any(part.height != cfg["expected_pool_rows"] for part in parts):
        raise ValueError("Incomplete model pool")
    if any(
        not parts[0][column].equals(parts[1][column])
        for column in ["compound_id", "fold", "y_true"]
    ):
        raise ValueError("Control and AP rows differ")
    candidate = parts[0].with_columns(
        pl.lit(cfg["output_model"]).alias("model"),
        ((1 - cfg["ap_weight"]) * parts[0]["score"] + cfg["ap_weight"] * parts[1]["score"]).alias(
            "score"
        ),
    )
    all_parts.extend([parts[0], candidate])
pl.concat(all_parts).write_parquet(outputs / "all-pool-oof.parquet")
