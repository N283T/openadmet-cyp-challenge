"""Refit the selected readout recipe on all endpoint labels and predict blinded test."""

import gc
import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(key):
        return pl.read_parquet(inp / cfg[key])

    labels = read("measurements").filter(
        pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
    )
    assignments = read("assignments")
    assert labels.height == 6525
    assert set(labels["compound_id"]).issubset(set(assignments["compound_id"]))
    aliases = (
        read("aliases")
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    ref = pl.read_csv(inp / cfg["reference"])
    test = (
        ref.select("SMILES", "Molecule_Name")
        .join(aliases, on="Molecule_Name", how="left", validate="1:1")
        .sort("compound_id")
    )
    assert test.height == test["compound_id"].n_unique() == 750
    assert test["compound_id"].null_count() == 0
    assert set(test["compound_id"]).isdisjoint(set(assignments["compound_id"]))
    emb = read("features")
    prefixes = cfg.get("prefixes", [cfg["prefix"]])
    cols = sorted(c for c in emb.columns if any(c.startswith(prefix) for prefix in prefixes))
    assert len(cols) == cfg["expected_dimensions"]
    assert emb["compound_id"].n_unique() == emb.height
    emb = emb.select("compound_id", *cols)
    auxcols = [f"aux_{i}" for i in range(4)] if cfg.get("auxiliary") else []
    if auxcols:
        aux = read("auxiliary")
        assert aux["compound_id"].n_unique() == aux.height
        emb = emb.join(aux.select("compound_id", *auxcols), on="compound_id", validate="1:1")
    assert (
        set(assignments["compound_id"]).union(test["compound_id"]).issubset(set(emb["compound_id"]))
    )
    threadpool_limits(limits=3)
    torch.set_num_threads(3)
    predictions, fits, members = [], [], []
    endpoints = cfg.get("endpoints", sorted(labels["endpoint"].unique()))
    assert endpoints and set(endpoints).issubset(set(labels["endpoint"].unique()))
    for endpoint in endpoints:
        train = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(emb, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        tt = test.join(emb, on="compound_id", validate="1:1", maintain_order="left")
        x, xt = (t.select(cols).to_numpy().astype(np.float32) for t in [train, tt])
        assert np.isfinite(x).all() and np.isfinite(xt).all()
        scaler = pca = None
        if cfg.get("pca"):
            scaler = StandardScaler().fit(x)
            pca = PCA(n_components=cfg["pca"], random_state=42, svd_solver="randomized").fit(
                scaler.transform(x)
            )
            x = pca.transform(scaler.transform(x)).astype(np.float32)
            xt = pca.transform(scaler.transform(xt)).astype(np.float32)
            assert int(scaler.n_samples_seen_) == train.height
        if auxcols:
            x = np.column_stack([x, train.select(auxcols).to_numpy().astype(np.float32)])
            xt = np.column_stack([xt, tt.select(auxcols).to_numpy().astype(np.float32)])
        assert np.isfinite(x).all() and np.isfinite(xt).all()
        y = train["value"].to_numpy()
        joblib.dump(
            dict(
                scaler=scaler,
                pca=pca,
                embedding_columns=cols,
                auxiliary_columns=auxcols,
                training_ids=train["compound_id"].to_list(),
            ),
            out / f"{endpoint}-preprocessing.joblib",
        )
        if cfg["backend"] == "tabpfn":
            from tabpfn import TabPFNRegressor
            from tabpfn.constants import ModelVersion

            model = TabPFNRegressor.create_default_for_version(
                ModelVersion(cfg.get("tabpfn_version", "v2.6")),
                n_estimators=cfg.get("tabpfn_estimators", 8),
                model_path=str(inp / cfg["weights"]),
                device="cuda",
                random_state=42,
                auto_scale_n_estimators=False,
                n_preprocessing_jobs=1,
                memory_saving_mode="auto",
                show_progress_bar=False,
                ignore_pretraining_limits=True,
                fit_mode="fit_preprocessors",
            )
            model.fit(x, y)
            pred = np.asarray(model.predict(xt), dtype=float)
            model.save_fit_state(out / f"{endpoint}.tabpfn_fit")
            model_count = 1
        elif cfg["backend"] == "tabicl":
            from tabicl import TabICLRegressor

            model = TabICLRegressor(
                n_estimators=8,
                batch_size=1,
                device="cuda",
                n_jobs=3,
                model_path=str(inp / cfg["weights"]),
                allow_auto_download=False,
                use_amp="auto",
                use_fa3="auto",
                offload_mode="auto",
                disk_offload_dir=os.environ["TMPDIR"],
                random_state=42,
            )
            model.fit(x, y)
            pred = np.asarray(model.predict(xt), dtype=float)
            joblib.dump(model, out / f"{endpoint}.tabicl.joblib")
            model_count = 1
        elif cfg["backend"] == "tuned_lgbm":
            summary = json.loads((inp / cfg["selection_summary"]).read_text())
            selections = sorted(
                [s for s in summary["selections"] if s["endpoint"] == endpoint],
                key=lambda s: s["fold"],
            )
            assert [s["fold"] for s in selections] == list(range(5))
            individual = []
            for s in selections:
                model = lgb.LGBMRegressor(
                    **s["params"],
                    n_jobs=3,
                    random_state=42 + s["fold"],
                    verbosity=-1,
                    deterministic=True,
                    force_col_wise=True,
                )
                model.fit(x, y)
                individual.append(model.predict(xt))
                model.booster_.save_model(str(out / f"{endpoint}-selection{s['fold']}.lgbm.txt"))
            pred = np.mean(individual, axis=0)
            model_count = 5
            (out / f"{endpoint}-selections.json").write_text(json.dumps(selections, indent=2))
        else:
            assert cfg["backend"] == "fixed_lgbm"
            model = lgb.LGBMRegressor(
                objective="regression_l1",
                n_estimators=300,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                reg_lambda=1.0,
                n_jobs=3,
                random_state=42,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            model.fit(x, y)
            pred = np.asarray(model.predict(xt), dtype=float)
            model.booster_.save_model(str(out / f"{endpoint}.lgbm.txt"))
            model_count = 1
        assert np.isfinite(pred).all()
        predictions.append(
            tt.select("compound_id", "Molecule_Name", "SMILES").with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.Series("prediction", pred),
                pl.lit("full_train").alias("split_id"),
                pl.lit(-1).alias("fold"),
                pl.lit(cfg["method"]).alias("method"),
            )
        )
        members.append(
            train.select("compound_id", "endpoint", "value").with_columns(
                pl.lit("full_training").alias("role"),
                pl.lit("full_train").alias("split_id"),
                pl.lit(-1).alias("fold"),
            )
        )
        fits.append(
            dict(
                endpoint=endpoint,
                train_n=train.height,
                test_n=tt.height,
                dimensions=x.shape[1],
                models=model_count,
                prediction_mean=float(pred.mean()),
                prediction_std=float(pred.std()),
            )
        )
        print(json.dumps(fits[-1]), flush=True)
        pl.concat(predictions).write_parquet(out / "test-predictions-progress.parquet")
        del model
        gc.collect()
        torch.cuda.empty_cache()
    pred_table = pl.concat(predictions)
    expected_test_rows = 750 * len(endpoints)
    assert (
        pred_table.height
        == pred_table.select("compound_id", "endpoint").n_unique()
        == expected_test_rows
    )
    pred_table.write_parquet(out / "test-predictions.parquet")
    membership = pl.concat(members)
    expected_membership = labels.filter(pl.col("endpoint").is_in(endpoints)).height
    assert membership.height == expected_membership
    membership.write_parquet(out / "membership.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                fits=fits,
                all_endpoint_labels_used=True,
                test_labels_used=False,
                calibration=False,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
