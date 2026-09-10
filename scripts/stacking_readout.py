"""TabPFN on train-only foundation PCA and nested out-of-sample log2fc predictions."""

import gc
import importlib.metadata
import json
import os
import time
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.decomposition import PCA
from split_baselines import metrics
from tabpfn import TabPFNClassifier, TabPFNRegressor
from tabpfn.constants import ModelVersion

from cyp_challenge.storage import digest, write_json


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)
    torch.manual_seed(config["seed"])
    fold = config["fold"]
    assignments = pl.read_parquet(next(inputs.glob("*-assignments.parquet"))).filter(
        pl.col("strategy") == "umap"
    )
    embedding = pl.concat(
        [pl.read_parquet(p) for p in sorted(inputs.glob("*-chemeleon-*.parquet"))]
    )
    feature_columns = sorted(c for c in embedding.columns if c.startswith("chemeleon_"))
    assert len(feature_columns) == 2048
    assert embedding["compound_id"].n_unique() == embedding.height == 6896
    auxiliary = []
    expected_heads = None
    for source in config["auxiliary"]:
        summary = json.loads((inputs / source["summary"]).read_text())
        assert summary["config"]["excluded_folds"] == source["excluded_folds"]
        assert fold in source["excluded_folds"]
        members = pl.read_parquet(inputs / source["membership"]).join(
            assignments.select("compound_id", "fold"), on="compound_id", validate="1:1"
        )
        assert members["fold"].null_count() == 0
        assert (
            members["fold"].is_in(source["excluded_folds"]) == (members["role"] == "outer_holdout")
        ).all()
        predictions = pl.concat([pl.read_parquet(inputs / p) for p in source["predictions"]])
        assert predictions["compound_id"].n_unique() == predictions.height == 6896
        heads = summary["heads"]
        assert len(heads) == 4
        if expected_heads is None:
            expected_heads = heads
        assert heads == expected_heads
        other = next(f for f in source["excluded_folds"] if f != fold)
        auxiliary.append((other, predictions.select("compound_id", *heads)))
    assert {a for a, _ in auxiliary} == set(range(5)) - {fold}
    labels = pl.read_parquet(next(inputs.glob("*-measurements.parquet")))
    aliases = pl.read_parquet(next(inputs.glob("*-aliases.parquet")))
    names = (
        aliases.select("compound_id", "molecule_name")
        .drop_nulls()
        .unique()
        .group_by("compound_id")
        .agg(pl.col("molecule_name").sort().str.join("|"))
    )
    checkpoints = {
        kind: next(inputs.glob(f"*-tabpfn-v2.6-{kind}-v2.6_default.ckpt"))
        for kind in ["regressor", "classifier"]
    }
    runtime = {
        "config": config,
        "versions": {
            n: importlib.metadata.version(n)
            for n in ["tabpfn", "torch", "scikit-learn", "numpy", "polars"]
        },
        "checkpoint_sha256": {k: digest(p) for k, p in checkpoints.items()},
        "cuda": torch.version.cuda,
        "device": torch.cuda.get_device_name(),
        "checkpoint_origin": "Existing local TabPFN-2.6 cache files, snapshotted and SHA256 pinned",
        "embedding_dimension": 260,
    }
    write_json(output / "runtime.json", runtime)
    all_predictions, scores = [], []
    for endpoint in config["endpoints"]:
        started = time.monotonic()
        frame = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(embedding, on="compound_id", how="left", validate="1:1")
            .join(
                assignments.select("compound_id", "fold"),
                on="compound_id",
                how="left",
                validate="1:1",
            )
            .join(names, on="compound_id", how="left", validate="1:1")
            .sort("compound_id")
        )
        assert frame["fold"].null_count() == 0
        x_raw = frame.select(feature_columns).to_numpy().astype(np.float32)
        assert np.isfinite(x_raw).all()
        y = frame["value"].to_numpy()
        train = frame["fold"].to_numpy() != fold
        val = ~train
        pca = PCA(n_components=config["pca_components"], random_state=config["seed"] + fold)
        reduced_train = pca.fit_transform(x_raw[train])
        reduced_val = pca.transform(x_raw[val])
        assert pca.n_samples_ == int(train.sum())
        np.testing.assert_allclose(pca.mean_, x_raw[train].mean(axis=0), atol=1e-6)
        write_json(
            output / f"{endpoint}-pca.json",
            {
                "fit_rows": int(train.sum()),
                "fit_fold_excluded": fold,
                "components": pca.n_components_,
                "explained_variance_ratio_sum": float(pca.explained_variance_ratio_.sum()),
                "heads": heads,
            },
        )
        joblib.dump(pca, output / f"{endpoint}-pca.joblib")
        aux = np.zeros((frame.height, 4), dtype=np.float32)
        for other, predictions in auxiliary:
            values = (
                frame.select("compound_id")
                .join(predictions, on="compound_id", how="left", validate="1:1")
                .select(heads)
                .to_numpy()
            )
            mask = frame["fold"].to_numpy() == other
            aux[mask] = values[mask]
            aux[val] += values[val] / 4
        x = np.empty((frame.height, 260), dtype=np.float32)
        x[train] = np.column_stack((reduced_train, aux[train]))
        x[val] = np.column_stack((reduced_val, aux[val]))
        assert np.isfinite(x).all()
        frame.select("compound_id", "fold").with_columns(
            pl.Series("features", x.tolist())
        ).write_parquet(output / f"{endpoint}-features.parquet")
        classification = endpoint.endswith("is_TDI")
        kind = "classifier" if classification else "regressor"
        model_type = TabPFNClassifier if classification else TabPFNRegressor
        model = model_type.create_default_for_version(
            ModelVersion.V2_6,
            model_path=str(checkpoints[kind]),
            device="cuda",
            n_estimators=config["n_estimators"],
            auto_scale_n_estimators=False,
            random_state=config["seed"] + fold,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
        )
        if classification:
            assert set(y[train]) == set(y[val]) == {0.0, 1.0}
        write_json(output / f"{endpoint}-parameters.json", model.get_params(deep=False))
        print(
            f"START {endpoint} fold={fold} train={int(train.sum())} "
            f"val={int(val.sum())} d={x.shape[1]}",
            flush=True,
        )
        model.fit(x[train], y[train])
        p = model.predict_proba(x[val])[:, 1] if classification else model.predict(x[val])
        assert np.isfinite(p).all()
        rows = (
            frame.filter(pl.Series(val))
            .select("compound_id", "molecule_name")
            .with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.lit(config["split_run"]).alias("split_id"),
                pl.lit("umap").alias("strategy"),
                pl.lit(fold).alias("fold"),
                pl.lit("nested_crossfit").alias("encoder_run"),
                pl.Series("y_true", y[val]),
                pl.Series("prediction", p),
            )
        )
        rows.write_parquet(output / f"{endpoint}-oof.parquet")
        all_predictions.append(rows)
        score = {
            "endpoint": endpoint,
            "fold": fold,
            "n": int(val.sum()),
            **metrics(
                y[val],
                p,
                frame["conf_low"].to_numpy()[val],
                frame["conf_high"].to_numpy()[val],
                classification,
            ),
        }
        scores.append(score)
        model.save_fit_state(output / f"{endpoint}.tabpfn_fit")
        print(json.dumps({**score, "seconds": time.monotonic() - started}), flush=True)
        del model, frame, x
        gc.collect()
        torch.cuda.empty_cache()
    pl.concat(all_predictions).write_parquet(output / "oof.parquet")
    pl.DataFrame(scores, infer_schema_length=None).write_parquet(output / "metrics.parquet")
    write_json(
        output / "summary.json",
        {**runtime, "metrics": scores, "oof_rows": sum(p.height for p in all_predictions)},
    )


if __name__ == "__main__":
    main()
