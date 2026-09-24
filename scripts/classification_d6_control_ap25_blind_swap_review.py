"""Review the exact four blind D6 swaps in the unsubmitted 75/25 probe."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
control = pl.read_parquet(inp / cfg["control"]).select(
    "Molecule_Name", pl.col(cfg["endpoint"]).alias("control_call")
)
candidate = pl.read_parquet(inp / cfg["candidate"]).select(
    "Molecule_Name", pl.col(cfg["endpoint"]).alias("candidate_call")
)
diagnostics = pl.read_parquet(inp / cfg["diagnostics"]).select(
    "molecule_name",
    "compound_id",
    "top10_similarity_weighted_positive_rate",
    "nearest_positive_similarity",
    "nearest_negative_similarity",
    "nearest_any_similarity",
)
decision = pl.read_parquet(inp / cfg["decision"]).select(
    "compound_id", pl.col("score").alias("blend_score")
)
review = (
    control.join(candidate, on="Molecule_Name", validate="1:1")
    .filter(pl.col("control_call") != pl.col("candidate_call"))
    .join(diagnostics, left_on="Molecule_Name", right_on="molecule_name", validate="1:1")
    .join(decision, on="compound_id", validate="1:1")
    .with_columns(
        pl.when(pl.col("candidate_call"))
        .then(pl.lit("added"))
        .otherwise(pl.lit("removed"))
        .alias("change")
    )
    .select(
        "Molecule_Name",
        "compound_id",
        "change",
        "blend_score",
        "top10_similarity_weighted_positive_rate",
        "nearest_positive_similarity",
        "nearest_negative_similarity",
        "nearest_any_similarity",
    )
    .sort("change", "Molecule_Name")
)
if review.height != 2 * cfg["expected_swaps_each_way"]:
    raise ValueError("Unexpected blind swap count")
review.write_csv(out / "blind-swaps.csv")
summary = review.group_by("change").agg(
    pl.len().alias("n"),
    pl.col("top10_similarity_weighted_positive_rate").mean().alias("mean_neighbor_positive_rate"),
    pl.col("nearest_positive_similarity").mean().alias("mean_nearest_positive_similarity"),
)
summary.write_csv(out / "swap-summary.csv")
print(summary, flush=True)
