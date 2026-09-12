"""Verify specialist OOF and memberships; evaluate fixed additions to global models."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    spec = importlib.util.spec_from_file_location(
        "official", next(inp.glob("*-custom_scoring_functions.py"))
    )
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    endpoint = "CYP2D6_pIC50_direct_inhibition"
    key = ["compound_id", "endpoint", "fold", "split_id"]
    global_oof = pl.read_parquet(inp / cfg["global_file"]).filter(pl.col("endpoint") == endpoint)
    assignments = pl.read_parquet(next(inp.glob("*-parent-assignments.parquet")))
    check = global_oof.join(
        assignments.select("compound_id", pl.col("fold").alias("assigned_fold"), "is_anchor"),
        on="compound_id",
        validate="m:1",
    )
    assert check.height == global_oof.height
    assert check.filter(
        (pl.col("fold") != pl.col("assigned_fold")) | pl.col("is_anchor")
    ).is_empty()
    # Global aggregate omits split_id; recover it only after verifying canonical folds.
    global_oof = global_oof.with_columns(pl.lit(cfg["split_run"]).alias("split_id"))
    baselines = {
        name: global_oof.filter(pl.col("method") == name).sort(key)
        for name in ["caruana_weight_crossfit", "caruana_full_fit"]
    }
    reference = baselines["caruana_weight_crossfit"]
    assert reference.height == 1457
    metrics, predictions, audits, correlations = [], [], [], []
    for name, item in cfg["models"].items():
        summary = json.loads((inp / item["summary"]).read_text())
        expected_heads = 16 if item["scope"] == "full" else 4 if item["scope"] == "d6_aux" else 1
        assert len(summary["heads"]) == expected_heads
        assert summary["config"]["selection_targets"] == ["CYP2D6"]
        for filename in item["memberships"]:
            m = pl.read_parquet(inp / filename)
            fold = m["outer_fold"].unique().item()
            outer = (pl.col("fold") == fold) & ~pl.col("is_anchor")
            assert m.filter(outer & (pl.col("role") != "outer_holdout")).is_empty()
            assert m.filter(~outer & (pl.col("role") == "outer_holdout")).is_empty()
            assert not set(m.filter(pl.col("role") == "training")["group"]) & set(
                m.filter(pl.col("role") == "inner_validation")["group"]
            )
            history = pl.read_ndjson(inp / item["histories"][str(fold)])
            head_column = history["validation_pic50_head_st_rae"].struct.field(endpoint)
            assert np.allclose(history["selection_score"], head_column, atol=1e-12, rtol=0)
            audits.append(
                dict(
                    model=name,
                    fold=fold,
                    heads=expected_heads,
                    exclusion_pass=True,
                    selection_target_verified=True,
                )
            )
        for backend in ["direct", "lgbm"]:
            table = (
                pl.read_parquet(inp / item[backend])
                .filter(pl.col("endpoint") == endpoint)
                .sort(key)
            )
            assert table.height == 1457 and table.select(key).equals(reference.select(key))
            for col in ["y_true", "conf_low", "conf_high"]:
                assert np.array_equal(table[col].to_numpy(), reference[col].to_numpy())
            assert np.isfinite(table["prediction"].to_numpy()).all()
            model = name + "_" + backend
            for baseline_name, baseline in baselines.items():
                p, b, y = (
                    table["prediction"].to_numpy(),
                    baseline["prediction"].to_numpy(),
                    table["y_true"].to_numpy(),
                )
                correlations.append(
                    dict(
                        model=model,
                        baseline=baseline_name,
                        prediction_pearson=float(pearsonr(p, b).statistic),
                        residual_pearson=float(pearsonr(p - y, b - y).statistic),
                        truth_spearman=float(spearmanr(p, y).statistic),
                    )
                )
                for weight in [0.0, 0.25, 0.5, 1.0]:
                    mixed = table.with_columns(
                        pl.Series("prediction", (1 - weight) * b + weight * p)
                    )
                    predictions.append(
                        mixed.select(
                            *key, "y_true", "conf_low", "conf_high", "prediction"
                        ).with_columns(
                            pl.lit(model).alias("model"),
                            pl.lit(baseline_name).alias("baseline"),
                            pl.lit(weight).alias("weight"),
                        )
                    )
                    for fold in [-1, 0, 1, 2, 3, 4]:
                        part = mixed if fold == -1 else mixed.filter(pl.col("fold") == fold)
                        score = scorer.rae_soft_threshold_absolute_error(
                            part["y_true"].to_numpy(),
                            part["prediction"].to_numpy(),
                            y_true_lower=part["conf_low"].to_numpy(),
                            y_true_upper=part["conf_high"].to_numpy(),
                        )
                        metrics.append(
                            dict(
                                model=model,
                                baseline=baseline_name,
                                weight=weight,
                                fold=fold,
                                n=part.height,
                                ST_RAE=float(score),
                            )
                        )
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    pl.DataFrame(correlations).write_csv(out / "correlations.csv")
    pl.concat(predictions).write_parquet(out / "fixed-blends.parquet")
    (out / "audit.json").write_text(
        json.dumps(dict(global_run=cfg["global_run"], checks=audits), indent=2)
    )
    print(
        pl.DataFrame(metrics).filter(
            (pl.col("fold") == -1) & (pl.col("baseline") == "caruana_weight_crossfit")
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
