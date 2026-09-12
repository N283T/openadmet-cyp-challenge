"""Audit TDI auxiliary readouts and their fixed additions to the current ensemble."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    spec = importlib.util.spec_from_file_location(
        "official", next(inp.glob("*-custom_scoring_functions.py"))
    )
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    keys = ["compound_id", "endpoint", "fold", "split_id"]
    baseline = pl.read_parquet(inp / cfg["baseline_file"]).sort(keys)
    rows = []
    predictions = []
    for method, filename in cfg["prediction_files"].items():
        table = pl.read_parquet(inp / filename).sort(keys)
        assert table.height == 6384 and table.select(keys).equals(baseline.select(keys))
        for col in ["y_true", "conf_low", "conf_high"]:
            assert np.array_equal(table[col].to_numpy(), baseline[col].to_numpy())
        for weight in [0.0, 0.25, 0.5, 1.0]:
            pred = (1 - weight) * baseline["prediction"].to_numpy() + weight * table[
                "prediction"
            ].to_numpy()
            mixed = table.with_columns(
                pl.Series("prediction", pred),
                pl.lit(method).alias("method"),
                pl.lit(weight).alias("candidate_weight"),
            )
            predictions.append(
                mixed.select(
                    *keys,
                    "y_true",
                    "conf_low",
                    "conf_high",
                    "is_neighbor",
                    "prediction",
                    "method",
                    "candidate_weight",
                )
            )
            for endpoint in sorted(table["endpoint"].unique()):
                part = mixed.filter(pl.col("endpoint") == endpoint)
                score = scorer.rae_soft_threshold_absolute_error(
                    part["y_true"].to_numpy(),
                    part["prediction"].to_numpy(),
                    y_true_lower=part["conf_low"].to_numpy(),
                    y_true_upper=part["conf_high"].to_numpy(),
                )
                rows.append(
                    dict(
                        method=method,
                        candidate_weight=weight,
                        endpoint=endpoint,
                        n=part.height,
                        st_rae=float(score),
                    )
                )
    metrics = pl.DataFrame(rows)
    metrics.write_csv(out / "metrics.csv")
    metrics.group_by("method", "candidate_weight").agg(
        pl.col("st_rae").mean().alias("macro_st_rae")
    ).sort("macro_st_rae").write_csv(out / "macro.csv")
    pl.concat(predictions).write_parquet(out / "prespecified-blends.parquet")
    audits = []
    for model, filenames in cfg["membership_files"].items():
        for filename in filenames:
            table = pl.read_parquet(inp / filename)
            fold = table["outer_fold"].unique().item()
            expected = (pl.col("fold") == fold) & ~pl.col("is_anchor")
            assert table.filter(expected & (pl.col("role") != "outer_holdout")).is_empty()
            assert table.filter(~expected & (pl.col("role") == "outer_holdout")).is_empty()
            assert table.filter(
                pl.col("is_anchor") & (pl.col("role") == "inner_validation")
            ).is_empty()
            train_groups = set(table.filter(pl.col("role") == "training")["group"])
            valid_groups = set(table.filter(pl.col("role") == "inner_validation")["group"])
            assert not train_groups & valid_groups
            audits.append(
                dict(
                    model=model,
                    fold=fold,
                    outer_compounds=table.filter(expected).height,
                    exclusion_passed=True,
                )
            )
    (out / "audit.json").write_text(
        json.dumps(dict(membership=audits, oof_rows=6384, baseline=cfg["baseline_run"]), indent=2)
    )
    print(metrics, flush=True)


if __name__ == "__main__":
    main()
