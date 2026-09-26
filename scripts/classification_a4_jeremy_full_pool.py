"""Extend the four saved Jeremy600 A4 readouts to the complete UMAP pool."""

import gc
import json
import os
import zipfile
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.preprocessing import FunctionTransformer
from tabicl import TabICLClassifier
from tabpfn.model_loading import load_fitted_tabpfn_model


def main() -> None:
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    scratch = Path(os.environ["TMPDIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)
    assignments = (
        pl.read_parquet(inputs / cfg["assignments"])
        .filter(pl.col("strategy") == cfg["split_id"])
        .select("compound_id", "fold", "group")
        .sort("compound_id")
    )
    labels = (
        pl.read_parquet(inputs / cfg["labels"])
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true"))
    )
    pool = assignments.join(labels, on="compound_id", how="left", validate="1:1")
    if pool.height != cfg["expected_pool_rows"]:
        raise ValueError("Incomplete UMAP pool")
    features = pl.read_parquet(inputs / cfg["features"])
    feature_columns = [column for column in features.columns if column.startswith("emb_")]
    raw = (
        pool.select("compound_id")
        .join(features, on="compound_id", validate="1:1")
        .select(feature_columns)
        .to_numpy()
    )
    if len(feature_columns) != 600 or not np.isfinite(raw).all():
        raise ValueError("Wrong Jeremy600 features")
    oof = {phase: pl.read_parquet(inputs / cfg[f"{phase}_oof"]) for phase in ("pca", "raw")}
    parts = []
    audits = []
    for view in cfg["views"]:
        phase, backend, method = view["phase"], view["backend"], view["method"]
        for fold in range(5):
            rows = np.flatnonzero(pool["fold"].to_numpy() == fold)
            held = pool.filter(pl.col("fold") == fold)
            prefix = f"{phase}-{backend}-fold{fold}"
            prep = joblib.load(inputs / f"{prefix}-prep.joblib")
            if prep["columns"] != feature_columns:
                raise ValueError("Feature columns changed")
            scaled = prep["scaler"].transform(raw[rows])
            x = (
                scaled.astype(np.float32)
                if prep["pca"] is None
                else prep["pca"].transform(scaled).astype(np.float32)
            )
            if backend == "tabicl":
                model = TabICLClassifier.load(inputs / f"{prefix}.fit", device="cuda")
            else:
                remapped = scratch / f"{prefix}.tabpfn_fit"
                with (
                    zipfile.ZipFile(inputs / f"{prefix}.fit") as source,
                    zipfile.ZipFile(remapped, "w") as dest,
                ):
                    for name in source.namelist():
                        payload = source.read(name)
                        if name == "init_params.json":
                            params = json.loads(payload)
                            params["model_path"] = str(inputs / cfg["pfn_weights"])
                            payload = json.dumps(params).encode()
                        dest.writestr(name, payload)
                model = load_fitted_tabpfn_model(remapped, device="cuda")
                if not hasattr(model, "date_transformer_"):
                    model.date_transformer_ = FunctionTransformer()
                if not hasattr(model, "text_transformer_"):
                    model.text_transformer_ = FunctionTransformer()
            probability = model.predict_proba(x)[:, list(model.classes_).index(1)]
            result = held.with_columns(
                pl.lit(cfg["endpoint"]).alias("endpoint"),
                pl.lit(cfg["split_id"]).alias("split_id"),
                pl.lit(view["name"]).alias("model"),
                pl.Series("score", probability.astype(float)),
            )
            labelled = result.filter(pl.col("y_true").is_not_null()).select("compound_id", "score")
            expected = oof[phase].filter(
                (pl.col("endpoint") == cfg["endpoint"])
                & (pl.col("method") == method)
                & (pl.col("fold") == fold)
            )
            aligned = labelled.join(
                expected.select("compound_id", pl.col("probability").alias("expected")),
                on="compound_id",
                validate="1:1",
            )
            if aligned.height != labelled.height:
                raise ValueError(f"Incomplete labelled replay for {prefix}")
            difference = np.abs(aligned["score"].to_numpy() - aligned["expected"].to_numpy())
            audits.append(
                {
                    "view": view["name"],
                    "fold": fold,
                    "labelled_rows": labelled.height,
                    "max_abs_diff": float(difference.max()),
                }
            )
            if difference.max() > 0.005:
                raise ValueError(f"OOF replay mismatch for {prefix}: {difference.max()}")
            parts.append(result)
            print(json.dumps(audits[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
    all_pool = pl.concat(parts).sort("model", "compound_id")
    if all_pool.height != pool.height * len(cfg["views"]):
        raise ValueError("Incomplete Jeremy all-pool bank")
    all_pool.write_parquet(outputs / "all-pool-oof.parquet")
    (outputs / "audit.json").write_text(json.dumps(audits, indent=2) + "\n")


if __name__ == "__main__":
    main()
