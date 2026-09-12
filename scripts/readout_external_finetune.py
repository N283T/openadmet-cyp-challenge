"""Fixed-CV readout with a different fine-tuned encoder for every outer fold.

feature_files maps string fold numbers to EXACT numbered snapshot filenames,
e.g. {"0": "0000-external-embeddings.parquet", ...}. Basenames may repeat.
"""

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
    paths = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
    assert len(paths) == 1, (name, paths)
    return pl.read_parquet(paths[0])


a, labels = read("parent-assignments.parquet"), read("measurements.parquet")
assert set(cfg["feature_files"]) == {str(fold) for fold in range(5)}
features = {}
feature_ids = None
for fold in range(5):
    filename = cfg["feature_files"][str(fold)]
    assert Path(filename).name == filename
    emb = pl.read_parquet(inp / filename)
    cols = sorted(c for c in emb.columns if c.startswith(cfg.get("feature_prefix", "emb_")))
    assert cols and emb["compound_id"].n_unique() == emb.height
    assert np.isfinite(emb.select(cols).to_numpy()).all()
    assert emb["fold"].unique().to_list() == [fold]
    assert emb["split_id"].unique().to_list() == [cfg["split_run"]]
    ids = set(emb["compound_id"])
    assert set(a["compound_id"]).issubset(ids)
    if feature_ids is not None:
        assert feature_ids == ids
    feature_ids = ids
    features[fold] = emb.select("compound_id", *cols)
torch.set_num_threads(4)
rows, fits = [], []
for ep in sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))["endpoint"].unique()
):
    base = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    for fold in range(5):
        emb = features[fold]
        cols = [c for c in emb.columns if c != "compound_id"]
        h = base.join(emb, on="compound_id", validate="1:1").sort("compound_id")
        assert h.height == base.height
        x, y = h.select(cols).to_numpy(), h["value"].to_numpy()
        va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
        xx = x
        if cfg.get("pca"):
            scaler = StandardScaler().fit(x[~va])
            pca = PCA(
                n_components=min(cfg["pca"], x.shape[1], int((~va).sum())), random_state=42
            ).fit(scaler.transform(x[~va]))
            xx = pca.transform(scaler.transform(x))
            joblib.dump((scaler, pca, cols), out / f"{ep}-{fold}-pca.joblib")
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
oof = pl.concat(rows)
expected = (
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
    .join(a, on="compound_id", validate="m:1")
    .filter(~pl.col("is_anchor"))
)
keys = ["compound_id", "endpoint", "fold"]
assert oof.height == expected.height == cfg.get("expected_oof_rows", 6384)
assert oof.select(keys).unique().height == oof.height
assert oof.select(keys).join(expected.select(keys), on=keys, how="anti").is_empty()
oof.write_parquet(out / "external-oof.parquet")
(out / "summary.json").write_text(json.dumps(dict(config=cfg, fits=fits), indent=2))
