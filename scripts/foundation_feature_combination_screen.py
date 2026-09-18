"""Evaluate selected heterogeneous feature combinations with foundation readouts."""

import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from tabicl import TabICLRegressor
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion
from threadpoolctl import threadpool_limits

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
torch.set_num_threads(cfg["num_threads"])
threadpool_limits(cfg["num_threads"])


def read(name: str) -> pl.DataFrame:
    path = inp / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pl.read_parquet(path)


spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
assignments = read("parent-assignments.parquet")
measurements = read("measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)

blocks = {}
block_columns = {}
for block_name, block_cfg in cfg["blocks"].items():
    frame = read(block_cfg["file"])
    source_columns = sorted(
        column
        for column in frame.columns
        if any(column.startswith(prefix) for prefix in block_cfg["prefixes"])
    )
    if not source_columns:
        raise ValueError(f"No columns for block {block_name}")
    renamed = {column: f"{block_name}::{column}" for column in source_columns}
    blocks[block_name] = frame.select("compound_id", *source_columns).rename(renamed)
    block_columns[block_name] = list(renamed.values())


def matrix(candidate: dict, table: pl.DataFrame):
    matrices = []
    names = []
    labels = []
    ids = table.select("compound_id")
    for block_name in candidate["blocks"]:
        aligned = ids.join(
            blocks[block_name],
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        columns = block_columns[block_name]
        values = aligned.select(columns).to_numpy().astype(np.float32)
        values[~np.isfinite(values)] = np.nan
        matrices.append(values)
        names.extend(columns)
        labels.extend([block_name] * len(columns))
    return np.column_stack(matrices), np.asarray(names), np.asarray(labels)


def select_features(x, y, training, names, fold, top_k):
    if not top_k:
        return np.arange(x.shape[1]), np.full(x.shape[1], np.nan)
    selector = LGBMRegressor(
        objective="regression_l1",
        n_estimators=300,
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.05,
        reg_lambda=1.0,
        n_jobs=cfg["num_threads"],
        random_state=cfg["seed"] + fold,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    selector.fit(x[training], y[training])
    gain = selector.booster_.feature_importance(importance_type="gain")
    order = np.lexsort((names, -gain))
    return order[: min(top_k, len(order))], gain


def finite_matrix(x, training):
    medians = np.nanmedian(x[training], axis=0)
    medians[~np.isfinite(medians)] = 0.0
    result = np.where(np.isfinite(x), x, medians).astype(np.float32)
    if not np.isfinite(result).all():
        raise ValueError("Non-finite values remain after training-median imputation")
    return result, medians


def make_tabpfn(fold: int):
    return TabPFNRegressor.create_default_for_version(
        ModelVersion(cfg["tabpfn_version"]),
        n_estimators=cfg["tabpfn_estimators"],
        model_path=str(inp / cfg["tabpfn_weights"]),
        device="cuda",
        random_state=cfg["seed"] + fold,
        auto_scale_n_estimators=False,
        n_preprocessing_jobs=1,
        memory_saving_mode="auto",
        show_progress_bar=False,
        ignore_pretraining_limits=True,
        fit_mode="fit_preprocessors",
    )


tabicl_params = {
    "n_estimators": cfg["tabicl_estimators"],
    "batch_size": 1,
    "device": "cuda",
    "n_jobs": cfg["num_threads"],
    "model_path": str(inp / cfg["tabicl_weights"]),
    "allow_auto_download": False,
    "use_amp": "auto",
    "use_fa3": "auto",
    "offload_mode": "auto",
    "disk_offload_dir": os.environ["TMPDIR"],
}
(out / "backend.json").write_text(
    json.dumps(
        {
            "tabpfn": {
                "version": importlib.metadata.version("tabpfn"),
                "model_version": cfg["tabpfn_version"],
                "estimators": cfg["tabpfn_estimators"],
            },
            "tabicl": {
                "version": importlib.metadata.version("tabicl"),
                "parameters": tabicl_params,
            },
        },
        indent=2,
    )
    + "\n"
)

oof_frames = []
fit_rows = []
metric_rows = []
selection_rows = []
for candidate in cfg["candidates"]:
    for endpoint in candidate["endpoints"]:
        table = (
            measurements.filter(pl.col("endpoint") == endpoint)
            .join(assignments, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x_raw, feature_names, feature_blocks = matrix(candidate, table)
        y = table["value"].to_numpy()
        folds = table["fold"].to_numpy()
        anchors = table["is_anchor"].to_numpy()
        for backend in candidate["backends"]:
            method = f"{candidate['name']}__{backend}"
            fold_frames = []
            for fold in range(5):
                validation = (folds == fold) & ~anchors
                training = ~validation
                selected, gain = select_features(
                    x_raw,
                    y,
                    training,
                    feature_names,
                    fold,
                    candidate.get("importance_top_k"),
                )
                x, medians = finite_matrix(x_raw[:, selected], training)
                if backend == "tabicl":
                    scaler = StandardScaler().fit(x[training])
                    x_model = scaler.transform(x).astype(np.float32)
                    model = TabICLRegressor(**tabicl_params, random_state=cfg["seed"] + fold)
                elif backend.startswith("tabpfn"):
                    scaler = None
                    x_model = x
                    model = make_tabpfn(fold)
                else:
                    raise ValueError(f"Unknown backend: {backend}")
                model.fit(x_model[training], y[training])
                prediction = np.asarray(model.predict(x_model[validation]), dtype=float)
                if not np.isfinite(prediction).all():
                    raise ValueError(f"Non-finite predictions for {method}/{endpoint}/{fold}")
                rows = (
                    table.filter(pl.Series(validation))
                    .select(
                        "compound_id",
                        "endpoint",
                        "fold",
                        "is_neighbor",
                        "conf_low",
                        "conf_high",
                        pl.col("value").alias("y_true"),
                    )
                    .with_columns(
                        pl.lit(method).alias("method"),
                        pl.lit(cfg["split_run"]).alias("split_id"),
                        pl.Series("prediction", prediction),
                    )
                )
                fold_frames.append(rows)
                fit_rows.append(
                    {
                        "method": method,
                        "endpoint": endpoint,
                        "fold": fold,
                        "input_dimensions": x_raw.shape[1],
                        "selected_dimensions": len(selected),
                        "training": int(training.sum()),
                        "validation": int(validation.sum()),
                        "imputation": "outer-training median",
                        "scaling": "outer-training standard"
                        if scaler is not None
                        else "backend internal",
                        "prediction_mean": float(prediction.mean()),
                    }
                )
                for rank, feature_index in enumerate(selected, start=1):
                    selection_rows.append(
                        {
                            "method": method,
                            "endpoint": endpoint,
                            "fold": fold,
                            "rank": rank,
                            "feature": feature_names[feature_index],
                            "block": feature_blocks[feature_index],
                            "selector_gain": None
                            if np.isnan(gain[feature_index])
                            else float(gain[feature_index]),
                            "training_median": float(medians[rank - 1]),
                        }
                    )
                print(json.dumps(fit_rows[-1]), flush=True)
                del model
                torch.cuda.empty_cache()
            endpoint_oof = pl.concat(fold_frames).sort("fold", "compound_id")
            endpoint_oof.write_parquet(out / f"{method}-{endpoint.split('_')[0]}-oof.parquet")
            oof_frames.append(endpoint_oof)
            prediction = endpoint_oof["prediction"].to_numpy()
            y_true = endpoint_oof["y_true"].to_numpy()
            metric_rows.append(
                {
                    "method": method,
                    "endpoint": endpoint,
                    "ST_RAE": float(
                        scorer.rae_soft_threshold_absolute_error(
                            y_true,
                            prediction,
                            y_true_lower=endpoint_oof["conf_low"].to_numpy(),
                            y_true_upper=endpoint_oof["conf_high"].to_numpy(),
                        )
                    ),
                    "MAE": float(np.abs(y_true - prediction).mean()),
                    "Spearman": float(spearmanr(y_true, prediction).statistic),
                    "n": endpoint_oof.height,
                }
            )

metrics = pl.DataFrame(metric_rows).sort("endpoint", "ST_RAE")
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(fit_rows).write_csv(out / "fits.csv")
pl.DataFrame(selection_rows).write_parquet(out / "selected-features.parquet")
pl.concat(oof_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps({"config": cfg, "metrics": metrics.to_dicts()}, indent=2) + "\n"
)
print(metrics, flush=True)
