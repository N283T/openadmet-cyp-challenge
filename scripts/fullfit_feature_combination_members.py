"""Full-fit selected heterogeneous feature members and predict blinded test."""

import gc
import importlib.metadata
import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
import torch
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


labels = read("measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
aliases = (
    read("aliases.parquet")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
    .sort("compound_id")
)
reference = pl.read_csv(inp / "reference.csv")
test = aliases.join(reference, on="Molecule_Name", validate="1:1", maintain_order="left")
if test.height != test["compound_id"].n_unique() or test.height != 750:
    raise ValueError("Invalid blinded-test identity table")

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


def matrix(candidate: dict, ids: pl.DataFrame):
    matrices = []
    names = []
    labels_ = []
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
        labels_.extend([block_name] * len(columns))
    return np.column_stack(matrices), np.asarray(names), np.asarray(labels_)


def fixed_lgbm(seed: int):
    return lgb.LGBMRegressor(
        objective="regression_l1",
        n_estimators=300,
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.05,
        reg_lambda=1.0,
        n_jobs=cfg["num_threads"],
        random_state=seed,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
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

prediction_frames = []
selection_rows = []
fit_rows = []
for candidate in cfg["candidates"]:
    endpoint = candidate["endpoint"]
    method = candidate["name"]
    train = labels.filter(pl.col("endpoint") == endpoint).sort("compound_id")
    train_ids = train.select("compound_id")
    test_ids = test.select("compound_id")
    x_raw, feature_names, feature_blocks = matrix(candidate, train_ids)
    xt_raw, test_feature_names, _ = matrix(candidate, test_ids)
    if not np.array_equal(feature_names, test_feature_names):
        raise ValueError(f"Train/test feature mismatch for {method}")
    y = train["value"].to_numpy()

    selected = np.arange(x_raw.shape[1])
    gain = np.full(x_raw.shape[1], np.nan)
    if candidate.get("importance_top_k"):
        selector = fixed_lgbm(cfg["seed"])
        selector.fit(x_raw, y)
        gain = selector.booster_.feature_importance(importance_type="gain")
        order = np.lexsort((feature_names, -gain))
        selected = order[: min(candidate["importance_top_k"], len(order))]
        selector.booster_.save_model(str(out / f"{method}-selector.lgbm.txt"))

    x = x_raw[:, selected]
    xt = xt_raw[:, selected]
    medians = np.nanmedian(x, axis=0)
    medians[~np.isfinite(medians)] = 0.0
    x_finite = np.where(np.isfinite(x), x, medians).astype(np.float32)
    xt_finite = np.where(np.isfinite(xt), xt, medians).astype(np.float32)
    scaler = None

    if candidate["backend"] == "lgbm":
        model = fixed_lgbm(cfg["seed"])
        model.fit(x, y)
        prediction = np.asarray(model.predict(xt), dtype=float)
        model.booster_.save_model(str(out / f"{method}.lgbm.txt"))
    elif candidate["backend"] == "tabpfn35fast":
        model = TabPFNRegressor.create_default_for_version(
            ModelVersion(cfg["tabpfn_version"]),
            n_estimators=cfg["tabpfn_estimators"],
            model_path=str(inp / cfg["tabpfn_weights"]),
            device="cuda",
            random_state=cfg["seed"],
            auto_scale_n_estimators=False,
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
        )
        model.fit(x_finite, y)
        prediction = np.asarray(model.predict(xt_finite), dtype=float)
        model.save_fit_state(out / f"{method}.tabpfn_fit")
    elif candidate["backend"] == "tabicl":
        scaler = StandardScaler().fit(x_finite)
        x_model = scaler.transform(x_finite).astype(np.float32)
        xt_model = scaler.transform(xt_finite).astype(np.float32)
        model = TabICLRegressor(**tabicl_params, random_state=cfg["seed"])
        model.fit(x_model, y)
        prediction = np.asarray(model.predict(xt_model), dtype=float)
        joblib.dump(model, out / f"{method}.tabicl.joblib")
    else:
        raise ValueError(f"Unknown backend: {candidate['backend']}")

    if not np.isfinite(prediction).all():
        raise ValueError(f"Non-finite prediction for {method}")
    joblib.dump(
        {
            "selected_indices": selected,
            "selected_features": feature_names[selected].tolist(),
            "training_medians": medians,
            "scaler": scaler,
            "training_ids": train["compound_id"].to_list(),
        },
        out / f"{method}-preprocessing.joblib",
    )
    for rank, feature_index in enumerate(selected, start=1):
        selection_rows.append(
            {
                "method": method,
                "endpoint": endpoint,
                "rank": rank,
                "feature": feature_names[feature_index],
                "block": feature_blocks[feature_index],
                "selector_gain": None
                if np.isnan(gain[feature_index])
                else float(gain[feature_index]),
                "training_median": float(medians[rank - 1]),
            }
        )
    prediction_frames.append(
        test.select("compound_id", "Molecule_Name", "SMILES").with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", prediction),
        )
    )
    fit_rows.append(
        {
            "method": method,
            "endpoint": endpoint,
            "backend": candidate["backend"],
            "training": train.height,
            "test": test.height,
            "input_dimensions": x_raw.shape[1],
            "selected_dimensions": len(selected),
            "prediction_mean": float(prediction.mean()),
            "prediction_std": float(prediction.std()),
        }
    )
    print(json.dumps(fit_rows[-1]), flush=True)
    del model
    gc.collect()
    torch.cuda.empty_cache()

predictions = pl.concat(prediction_frames).sort("endpoint", "method", "compound_id")
if predictions.select("compound_id", "endpoint", "method").n_unique() != predictions.height:
    raise ValueError("Duplicate full-fit prediction keys")
predictions.write_parquet(out / "test-predictions.parquet")
pl.DataFrame(selection_rows).write_parquet(out / "selected-features.parquet")
pl.DataFrame(fit_rows).write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "tabpfn_version": importlib.metadata.version("tabpfn"),
            "tabicl_version": importlib.metadata.version("tabicl"),
            "fits": fit_rows,
            "all_endpoint_labels_used": True,
            "test_labels_used": False,
        },
        indent=2,
    )
    + "\n"
)
