"""Recover canonical keyed predictions from immutable native all-head outputs."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / cfg["labels"])
    assignments = pl.read_parquet(inp / cfg["assignments"])
    parent = pl.read_parquet(inp / cfg["parent"])
    predictions = pl.read_parquet(inp / cfg["prediction_compounds"])
    mapping = pl.read_parquet(inp / cfg["test_mapping"])
    test = mapping.join(
        predictions.select("compound_id", "smiles"), on="compound_id", validate="1:1"
    )
    outputs, tests = [], []
    endpoints = (
        ["CYP3A4", "CYP2C9", "CYP2D6", "CYP1A2"]
        if cfg["recipe"] in ["aux8", "union18"]
        else ["CYP2D6"]
    )
    for fold in range(5):
        val = pl.read_parquet(inp / cfg[f"val{fold}"]).join(
            assignments.select("row_id", "compound_id", "effective_fold"),
            on="row_id",
            validate="1:1",
        )
        assert (val["effective_fold"] == fold).all()
        pred = pl.read_parquet(inp / cfg[f"test{fold}"]).join(
            test, on="compound_id", validate="1:1"
        )
        assert pred.height == 750
        for target in endpoints:
            head = target.lower() + "_pic50_direct_inhibition"
            endpoint = target + "_pIC50_direct_inhibition"
            native = val.select("compound_id", pl.col(head).alias("prediction"))
            scored = labels.filter(pl.col("endpoint") == endpoint).join(
                native, on="compound_id", validate="1:1"
            )
            expected = labels.filter(pl.col("endpoint") == endpoint).join(
                assignments.filter(pl.col("effective_fold") == fold).select("compound_id"),
                on="compound_id",
            )
            assert scored.height > 0 and scored.height == expected.height
            if cfg["cv_policy"] == "parent_adapted":
                canonical = (
                    labels.filter(pl.col("endpoint") == endpoint)
                    .join(parent, on="compound_id")
                    .filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
                )
                assert set(scored["compound_id"]) == set(canonical["compound_id"])
            outputs.append(
                scored.select(
                    "compound_id",
                    "endpoint",
                    pl.col("value").alias("y_true"),
                    "conf_low",
                    "conf_high",
                    "prediction",
                ).with_columns(
                    pl.lit(fold).alias("fold"),
                    pl.lit(cfg["split_id"]).alias("split_id"),
                    pl.lit(cfg["recipe"]).alias("method"),
                )
            )
            tests.append(
                pred.select(
                    "compound_id", "molecule_name", "smiles", pl.col(head).alias("prediction")
                ).with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(fold).alias("fold"),
                    pl.lit(cfg["recipe"]).alias("method"),
                )
            )
    oof = pl.concat(outputs)
    assert np.isfinite(oof["prediction"].to_numpy()).all()
    expected_rows = (
        (6384 if len(endpoints) == 4 else 1457) if cfg["cv_policy"] == "parent_adapted" else 1493
    )
    assert oof.height == expected_rows
    oof.write_parquet(out / "scp-oof.parquet")
    folds = pl.concat(tests)
    folds.write_parquet(out / "fold-test-predictions.parquet")
    mean = folds.group_by("compound_id", "molecule_name", "smiles", "endpoint", "method").agg(
        pl.col("prediction").mean(), pl.len().alias("n_folds")
    )
    assert mean["n_folds"].min() == mean["n_folds"].max() == 5
    mean.write_parquet(out / "test-predictions.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                source_training_run=cfg["source_run"],
                oof_rows=oof.height,
                test_rows=mean.height,
                fix="Correct canonical endpoint names when mapping saved native predictions; "
                "no model or prediction value changes",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
