"""Predict the complete UMAP pool with frozen TabPFN 2.6 fold fits."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
assigned = (
    pl.read_parquet(inp / cfg["assignments"])
    .filter(pl.col("strategy") == cfg["split_id"])
    .select("compound_id", "fold", "group")
)
labels = (
    pl.read_parquet(inp / cfg["labels"])
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true"))
)
features = pl.read_parquet(inp / cfg["features"])
if assigned.height != 6146 or assigned["compound_id"].n_unique() != 6146:
    raise ValueError("Incomplete assigned molecular pool")
pool = (
    assigned.join(labels, on="compound_id", how="left", validate="1:1")
    .join(features, on="compound_id", validate="1:1")
    .sort("compound_id")
)
torch.set_num_threads(2)
train = (
    labels.join(assigned, on="compound_id", validate="1:1")
    .join(features, on="compound_id", validate="1:1")
    .sort("compound_id")
)
parts = []
for fold in range(5):
    held = pool.filter(pl.col("fold") == fold)
    if "preprocess_aliases" in cfg:
        preprocess = joblib.load(inp / cfg["preprocess_aliases"][fold])
        columns = preprocess["columns"]
    else:
        columns = [
            column for column in features.columns if column.startswith(cfg["feature_prefix"])
        ]
        training = train.filter(pl.col("fold") != fold).select(columns).to_numpy()
        scaler = StandardScaler().fit(training)
        pca = (
            PCA(n_components=256, random_state=42).fit(scaler.transform(training))
            if cfg["preprocessing"] == "pca256"
            else None
        )
        preprocess = {"columns": columns, "scaler": scaler, "pca": pca}
    if len(columns) != 600 or not all(
        column.startswith(cfg["feature_prefix"]) for column in columns
    ):
        raise ValueError("Frozen feature columns differ")
    x = held.select(columns).to_numpy()
    if not np.isfinite(x).all():
        raise ValueError("Nonfinite feature")
    scaled = preprocess["scaler"].transform(x)
    x = (preprocess["pca"].transform(scaled) if preprocess["pca"] is not None else scaled).astype(
        np.float32
    )
    train_x = train.filter(pl.col("fold") != fold).select(columns).to_numpy()
    train_y = train.filter(pl.col("fold") != fold)["y_true"].to_numpy().astype(int)
    scaled_train = preprocess["scaler"].transform(train_x)
    train_x = (
        preprocess["pca"].transform(scaled_train) if preprocess["pca"] is not None else scaled_train
    ).astype(np.float32)
    model = TabPFNClassifier.create_default_for_version(
        ModelVersion.V2_6,
        n_estimators=cfg.get("n_estimators", 4),
        model_path=str(inp / cfg["weights"]),
        device="cuda",
        random_state=42 + fold,
        auto_scale_n_estimators=False,
        n_preprocessing_jobs=1,
        memory_saving_mode="auto",
        show_progress_bar=False,
        ignore_pretraining_limits=True,
        fit_mode="fit_preprocessors",
    )
    model.fit(train_x, train_y)
    positive = list(model.classes_).index(1)
    probability = model.predict_proba(x)[:, positive]
    if not np.isfinite(probability).all():
        raise ValueError(f"Nonfinite probability in fold {fold}")
    parts.append(
        held.select("compound_id", "fold", "group", "y_true").with_columns(
            pl.lit(cfg["endpoint"]).alias("endpoint"),
            pl.lit(cfg["split_id"]).alias("split_id"),
            pl.lit(cfg["model"]).alias("model"),
            pl.Series("score", probability),
        )
    )
    print(f"fold {fold} complete", flush=True)
result = pl.concat(parts).sort("compound_id")
if result.height != 6146 or result["compound_id"].n_unique() != 6146:
    raise ValueError("Incomplete prediction pool")
result.write_parquet(out / "all-pool-oof.parquet")
