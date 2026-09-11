"""TabICL compact/full baseline on exact saved parent-retained training features."""

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
from tabicl import TabICLRegressor
from threadpoolctl import threadpool_limits

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
fold = cfg["fold"]


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


a = read("parent-assignments")
labels = read("measurements")
emb = pl.concat([pl.read_parquet(p) for p in sorted(inp.glob("*-chemeleon-*.parquet"))])
desc = pl.concat([pl.read_parquet(p) for p in sorted(inp.glob("*-mordred-*.parquet"))])
ecols = sorted(c for c in emb.columns if c.startswith("chemeleon_"))
dcols = [c for c in desc.columns if c != "compound_id"]
auxiliary = []
for parent in cfg["auxiliary"]:
    meta = json.loads((inp / parent["summary"]).read_text())
    pair = meta["config"]["excluded_folds"]
    assert fold in pair
    members = pl.read_parquet(inp / parent["membership"]).join(
        a.select("compound_id", "fold"), on="compound_id", validate="1:1"
    )
    assert (members["fold"].is_in(pair) == (members["role"] == "outer_holdout")).all()
    pred = pl.concat([pl.read_parquet(inp / p) for p in parent["predictions"]])
    auxiliary.append((next(i for i in pair if i != fold), pred, meta["heads"]))
assert {g for g, _, _ in auxiliary} == set(range(5)) - {fold}
rows = []
audits = []
for ep in sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))["endpoint"].unique()
):
    f = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .join(desc, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    raw = f.select(ecols).to_numpy().astype(np.float32)
    dx = f.select(dcols).to_numpy().astype(float)
    dx[~np.isfinite(dx)] = np.nan
    y = f["value"].to_numpy()
    ff = f["fold"].to_numpy()
    anchor = f["is_anchor"].to_numpy()
    va = (ff == fold) & ~anchor
    aux = np.zeros((f.height, 4), dtype=np.float32)
    for other, pred, heads in auxiliary:
        values = (
            f.select("compound_id")
            .join(pred, on="compound_id", validate="1:1", maintain_order="left")
            .select(heads)
            .to_numpy()
        )
        aux[ff == other] = values[ff == other]
        aux[ff == fold] += values[ff == fold] / 4
    assert np.isfinite(aux).all()
    for mode in ["parent_retained"]:
        tr = (ff != fold) if mode == "parent_excluded" else ~va
        assert not (tr & va).any()
        if mode == "parent_retained":
            assert tr[anchor].all()
        saved = cfg["preprocessing"][ep]
        pca = joblib.load(inp / saved["pca"])
        imputer = joblib.load(inp / saved["imputer"])
        cols = json.loads((inp / saved["columns"]).read_text())
        with threadpool_limits(limits=4):
            reduced = pca.transform(raw)
        keep = [dcols.index(c) for c in cols]
        dd = imputer.transform(dx[:, keep])
        full = np.column_stack([reduced, aux, dd])
        booster = lgb.Booster(model_file=str(inp / saved["booster"]))
        reference = (
            read("parent-oof")
            .filter(
                (pl.col("endpoint") == ep)
                & (pl.col("mode") == mode)
                & (pl.col("features") == "full")
            )
            .sort("compound_id")
        )
        assert (
            reference["compound_id"].to_list() == f.filter(pl.Series(va))["compound_id"].to_list()
        )
        reconstructed = booster.predict(full[va], num_threads=4)
        if not np.allclose(reconstructed, reference["prediction"], atol=1e-10, rtol=1e-10):
            # Recreate the original fit/transform memory layout before checking again.
            with threadpool_limits(limits=4):
                restored = PCA(n_components=256, random_state=42 + fold).fit(raw[tr])
                rebuilt = restored.transform(raw)
            np.testing.assert_allclose(restored.components_, pca.components_, atol=1e-6)
            np.testing.assert_allclose(restored.mean_, pca.mean_, atol=1e-6)
            write_json(
                out / f"{ep}-pca-reconstruction.json",
                dict(
                    maximum_feature_difference=float(np.max(np.abs(rebuilt - reduced))),
                    maximum_initial_prediction_difference=float(
                        np.max(np.abs(reconstructed - reference["prediction"].to_numpy()))
                    ),
                ),
            )
            reduced = rebuilt
            full = np.column_stack([reduced, aux, dd])
        np.testing.assert_allclose(
            booster.predict(full[va], num_threads=4),
            reference["prediction"],
            atol=1e-10,
            rtol=1e-10,
        )
        xs = {"compact": np.column_stack([reduced, aux]), "full": full}
        for features, x in xs.items():
            torch.set_num_threads(4)
            model = TabICLRegressor(
                n_estimators=8,
                batch_size=1,
                device="cuda",
                random_state=42 + fold,
                n_jobs=4,
                model_path=str(inp / cfg["weights"]),
                allow_auto_download=False,
                use_amp="auto",
                use_fa3="auto",
                offload_mode="auto",
                disk_offload_dir=os.environ["TMPDIR"],
            )
            model.fit(x[tr], y[tr])
            p = model.predict(x[va])
            assert np.isfinite(p).all()

            rows.append(
                f.filter(pl.Series(va))
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
                    pl.lit(mode).alias("mode"),
                    pl.lit(features).alias("features"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", p),
                )
            )
            audits.append(
                dict(
                    endpoint=ep,
                    fold=fold,
                    mode=mode,
                    features=features,
                    train_n=int(tr.sum()),
                    val_n=int(va.sum()),
                    anchors_in_train=int((tr & anchor).sum()),
                    anchors_in_validation=int((va & anchor).sum()),
                    feature_count=x.shape[1],
                )
            )
            print(json.dumps(audits[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
pl.concat(rows).write_parquet(out / "tabicl-parent-oof.parquet")
write_json(out / "summary.json", dict(config=cfg, fit_audits=audits, aux_memberships_verified=True))
