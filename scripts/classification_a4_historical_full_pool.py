"""Extend the saved historical A4 TabPFN stack to the complete UMAP pool."""

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
from tabpfn.model_loading import load_fitted_tabpfn_model


def main() -> None:
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    scratch = Path(os.environ["TMPDIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    baseline = json.loads((inputs / cfg["baseline_config"]).read_text())
    def old(name: str) -> Path:
        return inputs / f"old-{name}"
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
    if pool.height != cfg["expected_pool_rows"] or labels.height != cfg["expected_labelled_rows"]:
        raise ValueError("Wrong molecular pool or label count")
    foundation = pl.concat([pl.read_parquet(old(name)) for name in baseline["foundation"]])
    descriptors = pl.concat([pl.read_parquet(old(name)) for name in baseline["mordred"]])
    foundation_columns = sorted(c for c in foundation.columns if c.startswith("chemeleon_"))
    raw = (
        pool.select("compound_id")
        .join(foundation, on="compound_id", validate="1:1")
        .select(foundation_columns)
        .to_numpy()
        .astype(np.float32)
    )
    if len(foundation_columns) != 2048 or not np.isfinite(raw).all():
        raise ValueError("Incomplete CheMeleon features")
    parts = []
    audits = []
    for item in baseline["folds"]:
        fold = item["fold"]
        held = pool.filter(pl.col("fold") == fold)
        rows = np.flatnonzero(pool["fold"].to_numpy() == fold)
        paths = item["endpoints"][cfg["endpoint"]]
        auxiliary = np.zeros((pool.height, 4), dtype=np.float32)
        for parent in item["auxiliary"]:
            meta = json.loads(old(parent["summary"]).read_text())
            if fold not in meta["config"]["excluded_folds"]:
                raise ValueError("Auxiliary encoder included the held fold")
            prediction = pl.concat([pl.read_parquet(old(name)) for name in parent["predictions"]])
            values = (
                pool.select("compound_id")
                .join(prediction, on="compound_id", validate="1:1")
                .select(meta["heads"])
                .to_numpy()
            )
            auxiliary += values / len(item["auxiliary"])
        columns = json.loads(old(paths["columns"]).read_text())
        raw_desc = (
            pool.select("compound_id")
            .join(descriptors, on="compound_id", validate="1:1")
            .select(columns)
            .to_numpy()
            .astype(np.float64)
        )
        raw_desc[~np.isfinite(raw_desc)] = np.nan
        imputer = joblib.load(old(paths["imputer"]))
        pca = joblib.load(old(paths["pca"]))
        x = np.column_stack(
            (pca.transform(raw[rows]), auxiliary[rows], imputer.transform(raw_desc[rows]))
        )
        if not np.isfinite(x).all():
            raise ValueError(f"Nonfinite model inputs in fold {fold}")
        remapped = scratch / f"a4-historical-fold{fold}.tabpfn_fit"
        with zipfile.ZipFile(old(paths["model"])) as source, zipfile.ZipFile(remapped, "w") as dest:
            for name in source.namelist():
                payload = source.read(name)
                if name == "init_params.json":
                    params = json.loads(payload)
                    params["model_path"] = str(old(baseline["weights"]["classifier"]))
                    payload = json.dumps(params).encode()
                dest.writestr(name, payload)
        model = load_fitted_tabpfn_model(remapped, device="cuda")
        # The current TabPFN predictor expects no-op date/text preprocessors that
        # were absent from this older numeric-only fit-state serialization.
        if not hasattr(model, "date_transformer_"):
            model.date_transformer_ = FunctionTransformer()
        if not hasattr(model, "text_transformer_"):
            model.text_transformer_ = FunctionTransformer()
        probability = model.predict_proba(x)[:, list(model.classes_).index(1)]
        result = held.with_columns(
            pl.lit(cfg["endpoint"]).alias("endpoint"),
            pl.lit(cfg["split_id"]).alias("split_id"),
            pl.lit(cfg["method"]).alias("model"),
            pl.Series("score", probability.astype(float)),
        )
        expected = pl.read_parquet(old(item["oof"])).filter(pl.col("endpoint") == cfg["endpoint"])
        labelled = result.filter(pl.col("y_true").is_not_null()).select("compound_id", "score")
        aligned = labelled.join(
            expected.select("compound_id", pl.col("prediction").alias("expected")),
            on="compound_id",
            validate="1:1",
        )
        if aligned.height != labelled.height:
            raise ValueError(f"OOF audit coverage mismatch in fold {fold}")
        difference = np.abs(aligned["score"].to_numpy() - aligned["expected"].to_numpy())
        audits.append(
            {
                "fold": fold,
                "labelled_rows": labelled.height,
                "max_abs_diff": float(difference.max()),
            }
        )
        if difference.max() > 0.005:
            raise ValueError(f"OOF replay mismatch in fold {fold}: {difference.max()}")
        parts.append(result)
        print(json.dumps(audits[-1]), flush=True)
        del model
        gc.collect()
        torch.cuda.empty_cache()
    full = pl.concat(parts).sort("compound_id")
    if full.height != pool.height or full["compound_id"].n_unique() != pool.height:
        raise ValueError("Incomplete historical all-pool OOF")
    full.write_parquet(outputs / "all-pool-oof.parquet")
    (outputs / "audit.json").write_text(json.dumps(audits, indent=2) + "\n")


if __name__ == "__main__":
    main()
