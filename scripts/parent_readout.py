"""Matched full-feature LightGBM readouts with parents retained versus excluded."""

import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
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
    for mode in ["parent_excluded", "parent_retained"]:
        tr = (ff != fold) if mode == "parent_excluded" else ~va
        assert not (tr & va).any()
        if mode == "parent_retained":
            assert tr[anchor].all()
        pca = PCA(n_components=256, random_state=42 + fold)
        with threadpool_limits(limits=4):
            pca.fit(raw[tr])
            reduced = pca.transform(raw)
        keep = np.isfinite(dx[tr]).any(axis=0)
        keep &= np.nan_to_num(np.nanstd(dx[tr], axis=0), nan=0) > 0
        imputer = SimpleImputer(strategy="median")
        imputer.fit(dx[tr][:, keep])
        dd = imputer.transform(dx[:, keep])
        write_json(
            out / f"{ep}-{mode}-columns.json", [c for c, k in zip(dcols, keep, strict=True) if k]
        )
        xs = {"mordred": dd, "full": np.column_stack([reduced, aux, dd])}
        joblib.dump(pca, out / f"{ep}-{mode}-pca.joblib")
        joblib.dump(imputer, out / f"{ep}-{mode}-imputer.joblib")
        for features, x in xs.items():
            model = lgb.LGBMRegressor(
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
            model.fit(x[tr], y[tr])
            p = model.predict(x[va])
            assert np.isfinite(p).all()
            model.booster_.save_model(str(out / f"{ep}-{mode}-{features}.txt"))
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
pl.concat(rows).write_parquet(out / "parent-oof.parquet")
write_json(out / "summary.json", dict(config=cfg, fit_audits=audits, aux_memberships_verified=True))
