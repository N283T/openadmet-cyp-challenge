"""Target-specific fixed-split readout screen over multiple checkpoint views."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from catboost import CatBoostRegressor
from scipy.stats import spearmanr
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
torch.set_num_threads(cfg["num_threads"])

assignments = pl.read_parquet(inp / "parent-assignments.parquet")
measurements = pl.read_parquet(inp / "measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def make_model(backend: str, fold: int):
    if backend == "lgbm":
        return lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.05,
            reg_lambda=1.0,
            n_jobs=cfg["num_threads"],
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
    if backend == "catboost":
        return CatBoostRegressor(
            iterations=500,
            depth=6,
            learning_rate=0.03,
            l2_leaf_reg=5,
            loss_function="MAE",
            random_seed=42 + fold,
            thread_count=cfg["num_threads"],
            verbose=False,
            allow_writing_files=False,
            task_type="CPU",
        )
    if backend == "tabpfn":
        version = ModelVersion(cfg.get("tabpfn_version", "v2.6"))
        return TabPFNRegressor.create_default_for_version(
            version,
            n_estimators=cfg.get("tabpfn_estimators", 8),
            model_path=str(inp / cfg["tabpfn_weights"]),
            device="cuda",
            random_state=42 + fold,
            auto_scale_n_estimators=False,
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
        )
    raise ValueError(f"Unknown backend: {backend}")


all_oof = []
metric_rows = []
fit_rows = []
for candidate in cfg["candidates"]:
    features = pl.read_parquet(inp / candidate["file"])
    columns = sorted(
        column
        for column in features.columns
        if any(column.startswith(prefix) for prefix in candidate["prefixes"])
    )
    assert columns and features["compound_id"].n_unique() == features.height
    assert np.isfinite(features.select(columns).to_numpy()).all()
    nested = (
        pl.read_parquet(inp / candidate["nested_file"])
        if candidate.get("nested_file")
        else None
    )
    for backend in candidate["backends"]:
        method = f"{candidate['name']}__{backend}"
        endpoint_oof = []
        for endpoint in sorted(measurements["endpoint"].unique()):
            table = (
                measurements.filter(pl.col("endpoint") == endpoint)
                .join(assignments, on="compound_id", validate="1:1")
                .join(features, on="compound_id", validate="1:1")
                .sort("compound_id")
            )
            x = table.select(columns).to_numpy()
            y = table["value"].to_numpy()
            fold_rows = []
            for fold in range(5):
                validation = (table["fold"].to_numpy() == fold) & ~table[
                    "is_anchor"
                ].to_numpy()
                xx = x
                if nested is not None:
                    auxiliary_columns = [f"aux_{index}" for index in range(4)]
                    auxiliary = table.select("compound_id").join(
                        nested.filter(pl.col("outer_fold") == fold),
                        on="compound_id",
                        how="left",
                        validate="1:1",
                        maintain_order="left",
                    )
                    assert auxiliary.select(auxiliary_columns).null_count().sum_horizontal()[0] == 0
                    xx = np.column_stack(
                        [x, auxiliary.select(auxiliary_columns).to_numpy()]
                    )
                model = make_model(backend, fold)
                model.fit(xx[~validation], y[~validation])
                prediction = np.asarray(model.predict(xx[validation]), dtype=float)
                assert np.isfinite(prediction).all()
                rows = table.filter(pl.Series(validation)).select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                ).with_columns(
                    pl.lit(method).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", prediction),
                )
                fold_rows.append(rows)
                fit_rows.append(
                    {
                        "method": method,
                        "endpoint": endpoint,
                        "fold": fold,
                        "dimensions": xx.shape[1],
                        "train": int((~validation).sum()),
                        "validation": int(validation.sum()),
                    }
                )
                print(json.dumps(fit_rows[-1]), flush=True)
                del model
                if backend == "tabpfn":
                    torch.cuda.empty_cache()
            endpoint_frame = pl.concat(fold_rows)
            endpoint_oof.append(endpoint_frame)
            y_true = endpoint_frame["y_true"].to_numpy()
            prediction = endpoint_frame["prediction"].to_numpy()
            metric_rows.append(
                {
                    "method": method,
                    "endpoint": endpoint,
                    "ST_RAE": float(
                        scorer.rae_soft_threshold_absolute_error(
                            y_true,
                            prediction,
                            y_true_lower=endpoint_frame["conf_low"].to_numpy(),
                            y_true_upper=endpoint_frame["conf_high"].to_numpy(),
                        )
                    ),
                    "MAE": float(np.abs(y_true - prediction).mean()),
                    "Spearman": float(spearmanr(y_true, prediction).statistic),
                    "n": endpoint_frame.height,
                }
            )
        candidate_oof = pl.concat(endpoint_oof).sort("endpoint", "fold", "compound_id")
        candidate_oof.write_parquet(out / f"{method}-oof.parquet")
        all_oof.append(candidate_oof)

metrics = pl.DataFrame(metric_rows).sort("endpoint", "ST_RAE")
macro = metrics.group_by("method").agg(
    pl.col("ST_RAE", "MAE", "Spearman").mean()
).sort("ST_RAE")
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
pl.DataFrame(fit_rows).write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "primary_selection_scope": "endpoint",
            "macro_use": "descriptive_only",
            "endpoint_metrics": metrics.to_dicts(),
            "macro": macro.to_dicts(),
        },
        indent=2,
    )
)
print(macro)
print(metrics)
