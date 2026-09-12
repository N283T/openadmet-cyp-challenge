"""Canonical classification labels and fixed original UMAP membership, no model fitting."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
paths = [next(inp.glob(f"{i:03d}-*")) for i in range(5)]
measurements, assignments, compounds, test = [pl.read_parquet(p) for p in paths[:4]]
assignments = assignments.filter(pl.col("strategy") == "umap")
assert assignments["compound_id"].n_unique() == assignments.height == 6146
assert set(assignments["fold"]) == set(range(5))
assert set(assignments["compound_id"]).isdisjoint(set(test["compound_id"]))
labels = measurements.filter(pl.col("endpoint").is_in(["CYP2D6_is_TDI", "CYP3A4_is_TDI"]))
assert labels.height == 5081 and set(labels["value"]) == {0.0, 1.0}
assert not labels.select("compound_id", "endpoint").is_duplicated().any()
labels = labels.join(assignments, on="compound_id", validate="m:1").with_columns(
    pl.col("value").cast(pl.Int8).alias("y_true"),
    pl.lit("20260909T135130268493Z-edba6fc0").alias("split_id"),
)
assert labels.height == 5081 and labels["fold"].null_count() == 0
labels.write_parquet(out / "classification-labels.parquet")
assignments.with_columns(pl.lit("20260909T135130268493Z-edba6fc0").alias("split_id")).write_parquet(
    out / "umap-assignments.parquet"
)
counts = (
    labels.group_by("endpoint", "fold")
    .agg(
        pl.len().alias("n"),
        pl.col("y_true").sum().alias("positives"),
        pl.col("y_true").mean().alias("positive_rate"),
    )
    .sort("endpoint", "fold")
)
counts.with_columns((pl.col("n") - pl.col("positives")).alias("negatives")).write_csv(
    out / "fold-balance.csv"
)
overall = (
    labels.group_by("endpoint")
    .agg(
        pl.len().alias("n"),
        pl.col("y_true").sum().alias("positives"),
        pl.col("y_true").mean().alias("positive_rate"),
    )
    .sort("endpoint")
)
overall.write_csv(out / "label-balance.csv")
joined = compounds.join(
    labels.select("compound_id", "endpoint", "y_true").pivot(
        on="endpoint", index="compound_id", values="y_true"
    ),
    on="compound_id",
    how="left",
    validate="1:1",
)
joined.write_parquet(out / "classification-compounds.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            label_rows=labels.height,
            compound_rows=compounds.height,
            assigned_train_compounds=assignments.height,
            test_compounds=test.height,
            classification_labels="Use supplied boolean labels without rederivation; "
            "ambiguity remains missing",
            metric="MCC on confident CYP2D6_is_TDI and CYP3A4_is_TDI calls",
            source_policy="Use frozen external embeddings/static descriptors; "
            "no observed endpoint/TDI/pIC50/shift inputs. Challenge-finetuned predictions "
            "require classification-UMAP exclusion of all heads.",
            old_regression_parent="a3923255 is NOT this classification UMAP split",
            label_balance=overall.to_dicts(),
            assignment_columns=assignments.columns,
        ),
        indent=2,
    )
)
print(overall.write_csv())
print(counts.write_csv())
