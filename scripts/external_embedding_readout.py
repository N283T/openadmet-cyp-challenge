"""Fixed-CV external representation readout; all inputs explicitly snapshotted."""

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
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    exact = inp / name
    paths = [exact] if exact.exists() else [
        p for p in inp.iterdir() if p.name.partition("-")[2] == name
    ]
    assert len(paths) == 1, (name, paths)
    return pl.read_parquet(paths[0])


a, labels = read("parent-assignments.parquet"), read("measurements.parquet")
emb = read(cfg["feature_file"])
nested = read(cfg["nested_file"]) if cfg.get("nested_file") else None
cols = sorted(c for c in emb.columns if c.startswith(cfg.get("feature_prefix", "emb_")))
assert cols and emb["compound_id"].n_unique() == emb.height
assert np.isfinite(emb.select(cols).to_numpy()).all()
torch.set_num_threads(4)
rows, fits = [], []
available_endpoints = sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))["endpoint"].unique()
)
endpoints = cfg.get("endpoints", available_endpoints)
assert endpoints and set(endpoints).issubset(set(available_endpoints))
for ep in endpoints:
    h = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    x, y = h.select(cols).to_numpy(), h["value"].to_numpy()
    for fold in cfg.get("folds", list(range(5))):
        va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
        xx = x
        if cfg.get("mask_observed_validation"):
            xx = x.copy()
            for j, name in enumerate(cols):
                if name.startswith("score_fused_"):
                    original = name.replace("score_fused_", "score_", 1)
                    xx[va, j] = x[va, cols.index(original)]
                elif name.startswith("score_available_"):
                    xx[va, j] = 0
        if cfg.get("pca"):
            scaler = StandardScaler().fit(xx[~va])
            pca = PCA(
                n_components=min(cfg["pca"], x.shape[1], int((~va).sum())), random_state=42
            ).fit(scaler.transform(xx[~va]))
            xx = pca.transform(scaler.transform(xx))
            joblib.dump((scaler, pca, cols), out / f"{ep}-{fold}-pca.joblib")
        if nested is not None:
            auxiliary = h.select("compound_id").join(
                nested.filter(pl.col("outer_fold") == fold),
                on="compound_id",
                how="left",
                validate="1:1",
                maintain_order="left",
            )
            acols = [f"aux_{i}" for i in range(4)]
            assert all(c in auxiliary.columns for c in acols)
            ax = auxiliary.select(acols).to_numpy()
            assert np.isfinite(ax).all()
            xx = np.column_stack([xx, ax])
        if cfg["backend"] == "lgbm":
            params = dict(
                objective="regression_l1",
                n_estimators=300,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                reg_lambda=1,
                n_jobs=4,
                random_state=42 + fold,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            params.update(cfg.get("lgbm_overrides", {}))
            model = lgb.LGBMRegressor(**params)
        else:
            assert cfg["backend"] == "tabpfn"
            model = TabPFNRegressor.create_default_for_version(
                ModelVersion.V2_6,
                n_estimators=8,
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
        model.fit(xx[~va], y[~va])
        pred = model.predict(xx[va]).astype(float)
        assert np.isfinite(pred).all()
        if cfg["backend"] == "lgbm":
            model.booster_.save_model(str(out / f"{ep}-{fold}.lgbm.txt"))
        else:
            model.save_fit_state(out / f"{ep}-{fold}.tabpfn_fit")
        rows.append(
            h.filter(pl.Series(va))
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
                pl.lit(cfg["method"]).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.Series("prediction", pred),
            )
        )
        fits.append(
            dict(
                endpoint=ep,
                fold=fold,
                train=int((~va).sum()),
                validation=int(va.sum()),
                dimensions=xx.shape[1],
            )
        )
        print(json.dumps(fits[-1]), flush=True)
pl.concat(rows).write_parquet(out / "external-oof.parquet")
(out / "summary.json").write_text(json.dumps(dict(config=cfg, fits=fits), indent=2))
