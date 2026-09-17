"""TabICL on frozen or fold-specific external embeddings with optional nested log2fc.

Fit standardization and optional PCA to embedding columns using outer training
rows only; append the four nested auxiliary columns after embedding preprocessing.
Set ``pca`` to zero to retain every standardized embedding coordinate.
"""

import importlib.metadata
import inspect
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tabicl import TabICLRegressor
from threadpoolctl import threadpool_limits


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(inp / name)

    a, labels = read(cfg["assignments"]), read(cfg["measurements"])
    nested = read(cfg["nested_file"]) if cfg.get("nested_file") else None
    features = {}
    feature_ids = None
    prefixes = cfg.get("feature_prefix", "emb_")
    assert set(cfg["feature_files"]) == set(map(str, range(5)))
    for fold in range(5):
        emb = read(cfg["feature_files"][str(fold)])
        cols = sorted(c for c in emb.columns if c.startswith(prefixes))
        assert cols and not set(cols) & set(["value", "conf_low", "conf_high", "prediction"])
        assert emb["compound_id"].n_unique() == emb.height
        if cfg["fold_specific"]:
            assert emb["fold"].unique().to_list() == [fold]
            assert emb["split_id"].unique().to_list() == [cfg["split_run"]]
        ids = set(emb["compound_id"])
        assert set(a["compound_id"]).issubset(ids)
        if feature_ids is not None:
            assert ids == feature_ids
        feature_ids = ids
        assert np.isfinite(emb.select(cols).to_numpy()).all()
        features[fold] = emb.select("compound_id", *cols)
    auxcols = [f"aux_{i}" for i in range(4)]
    if nested is not None:
        assert set(nested["outer_fold"].unique()) == set(range(5))
        assert nested.select("compound_id", "outer_fold").n_unique() == nested.height
        for fold in range(5):
            aux = nested.filter(pl.col("outer_fold") == fold)
            assert set(a["compound_id"]).issubset(set(aux["compound_id"]))
            assert np.isfinite(aux.select(auxcols).to_numpy()).all()
    threadpool_limits(limits=cfg.get("threads", 2))
    torch.set_num_threads(cfg.get("threads", 2))
    common_params = dict(
        n_estimators=cfg.get("n_estimators", 8),
        batch_size=1,
        device=cfg.get("device", "cuda"),
        n_jobs=cfg.get("threads", 2),
        model_path=str(inp / cfg["weights"]),
        allow_auto_download=False,
        use_amp="auto",
        use_fa3="auto",
        offload_mode="auto",
        disk_offload_dir=os.environ["TMPDIR"],
    )
    # Real installed constructor validation does not fit or initialize a GPU model.
    constructor = TabICLRegressor(**common_params, random_state=42)
    assert constructor.get_params()["allow_auto_download"] is False
    (out / "backend.json").write_text(
        json.dumps(
            dict(
                version=importlib.metadata.version("tabicl"),
                signature=str(inspect.signature(TabICLRegressor)),
                parameters=common_params,
            ),
            indent=2,
        )
    )
    rows, fits, members = [], [], []
    for (endpoint,), base in labels.filter(
        pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
    ).group_by("endpoint", maintain_order=True):
        base = base.join(a, on="compound_id", validate="1:1").sort("compound_id")
        for fold in range(5):
            emb = features[fold]
            cols = [c for c in emb.columns if c != "compound_id"]
            h = base.join(emb, on="compound_id", validate="1:1", maintain_order="left")
            assert h.height == base.height and h["compound_id"].equals(base["compound_id"])
            va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
            assert va.any() and (~va).any()
            raw = h.select(cols).to_numpy().astype(np.float32)
            scaler = StandardScaler().fit(raw[~va])
            requested_pca = cfg.get("pca", 256)
            assert requested_pca >= 0
            scaled = scaler.transform(raw)
            if requested_pca:
                ncomp = min(requested_pca, raw.shape[1], int((~va).sum()) - 1)
                assert ncomp > 0
                pca = PCA(n_components=ncomp, random_state=42, svd_solver="randomized").fit(
                    scaled[~va]
                )
                x = pca.transform(scaled).astype(np.float32)
            else:
                ncomp = raw.shape[1]
                pca = None
                x = scaled.astype(np.float32)
            assert x.shape[1] == ncomp
            if nested is not None:
                aux = h.select("compound_id").join(
                    nested.filter(pl.col("outer_fold") == fold).select("compound_id", *auxcols),
                    on="compound_id",
                    validate="1:1",
                    maintain_order="left",
                )
                assert aux.height == h.height and aux["compound_id"].equals(h["compound_id"])
                values = aux.select(auxcols).to_numpy().astype(np.float32)
                x = np.column_stack([x, values])
                assert np.array_equal(x[:, -4:], values)
            assert np.isfinite(x).all()
            joblib.dump(
                dict(
                    scaler=scaler,
                    pca=pca,
                    embedding_columns=cols,
                    auxiliary_columns=auxcols if nested is not None else [],
                    preprocessing_training_compound_ids=h.filter(pl.Series(~va))[
                        "compound_id"
                    ].to_list(),
                ),
                out / f"{endpoint}-fold{fold}-preprocessing.joblib",
            )
            members.append(
                h.select("compound_id", "endpoint", "group", "is_anchor").with_columns(
                    pl.lit(fold).alias("outer_fold"),
                    pl.Series("role", np.where(va, "outer_validation", "outer_training")),
                )
            )
            meta = dict(
                endpoint=endpoint,
                fold=fold,
                training_n=int((~va).sum()),
                validation_n=int(va.sum()),
                embedding_columns=len(cols),
                pca_components=0 if pca is None else ncomp,
                embedding_output_columns=ncomp,
                auxiliary_columns=4 if nested is not None else 0,
                final_columns=x.shape[1],
                pca_n_samples=int(scaler.n_samples_seen_),
                feature_snapshot=cfg["feature_files"][str(fold)],
                nested_snapshot=cfg.get("nested_file"),
                embedding_preprocessing="standard_scaler_only"
                if pca is None
                else "standard_scaler_then_pca",
                method=cfg["method"],
            )
            assert meta["pca_n_samples"] == meta["training_n"]
            if cfg.get("prepare_only", False):
                print(json.dumps(meta), flush=True)
            else:
                model = TabICLRegressor(**common_params, random_state=42 + fold)
                model.fit(x[~va], h["value"].to_numpy()[~va])
                prediction = np.asarray(model.predict(x[va]), dtype=float)
                assert np.isfinite(prediction).all()
                joblib.dump(model, out / f"{endpoint}-fold{fold}.tabicl.joblib")
                rows.append(
                    h.filter(pl.Series(va))
                    .select(
                        "compound_id",
                        "endpoint",
                        "fold",
                        "is_neighbor",
                        pl.col("value").alias("y_true"),
                        "conf_low",
                        "conf_high",
                    )
                    .with_columns(
                        pl.Series("prediction", prediction),
                        pl.lit(cfg["method"]).alias("method"),
                        pl.lit(cfg["split_run"]).alias("split_id"),
                    )
                )
                print(json.dumps(meta | dict(prediction_mean=float(prediction.mean()))), flush=True)
            fits.append(meta)
    assert sum(r["validation_n"] for r in fits) == 6384
    pl.concat(members).write_parquet(out / "readout-membership.parquet")
    (out / "fits.json").write_text(json.dumps(fits, indent=2))
    if not cfg.get("prepare_only", False):
        oof = pl.concat(rows).sort("compound_id", "endpoint", "fold")
        assert (
            oof.height == 6384 and oof.select("compound_id", "endpoint", "fold").n_unique() == 6384
        )
        oof.write_parquet(out / "external-oof.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                prepared_fits=len(fits),
                oof_rows=0 if cfg.get("prepare_only", False) else 6384,
                expected_oof_rows=6384,
                fold_specific=cfg["fold_specific"],
                features=(
                    "Train-only embedding scaling with optional PCA, then four optional keyed "
                    "nested log2fc values"
                ),
                status="prepared_only" if cfg.get("prepare_only", False) else "completed",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
