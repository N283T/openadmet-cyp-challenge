"""Reconstruct honest inner predictions for the historical nested TabPFN stack."""

import gc
import json
import os
import time
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
assign = pl.read_parquet(inp / cfg["assignments"]).filter(pl.col("strategy") == "umap")
labels = pl.read_parquet(inp / cfg["labels"])
emb = pl.concat([pl.read_parquet(inp / name) for name in cfg["embedding"]])
ecols = sorted(c for c in emb.columns if c.startswith("chemeleon_"))
desc = pl.concat([pl.read_parquet(inp / name) for name in cfg["mordred"]])
dcols = [c for c in desc.columns if c != "compound_id"]
assert len(ecols) == 2048
endpoints = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
old_outer = pl.read_parquet(inp / cfg["old_outer"]).filter(pl.col("endpoint").is_in(endpoints))
old_test = pl.read_parquet(inp / cfg["old_test"]).filter(pl.col("endpoint").is_in(endpoints))
assert old_outer.height == 5081 and old_test.height == 7500
canonical = labels.filter(pl.col("endpoint").is_in(endpoints)).join(
    assign.select("compound_id", "fold"), on="compound_id", validate="m:1"
)
paired = old_outer.join(
    canonical.select("compound_id", "endpoint", "fold", "value"),
    on=["compound_id", "endpoint", "fold"],
    validate="1:1",
)
assert paired.height == 5081
np.testing.assert_array_equal(paired["y_true"], paired["value"])
old_outer.select(
    "compound_id", "endpoint", "fold", "y_true", pl.col("prediction").alias("probability")
).with_columns(
    pl.lit("historical_stack_mordred").alias("method"), pl.lit(cfg["split_id"]).alias("split_id")
).write_parquet(out / "outer-oof.parquet")
old_test.select(
    "compound_id", "endpoint", "fold", pl.col("prediction").alias("probability")
).with_columns(
    pl.lit("historical_stack_mordred").alias("method"), pl.lit(cfg["split_id"]).alias("split_id")
).write_parquet(out / "fold-test-probabilities.parquet")
old_mean = pl.read_parquet(inp / cfg["old_mean"]).filter(pl.col("endpoint").is_in(endpoints))
rebuilt = old_test.group_by("compound_id", "endpoint").agg(
    pl.col("prediction").mean().alias("rebuilt"), pl.len().alias("count")
)
paired_mean = old_mean.join(rebuilt, on=["compound_id", "endpoint"], validate="1:1")
assert paired_mean.height == 1500 and (paired_mean["count"] == 5).all()
np.testing.assert_allclose(paired_mean["prediction"], paired_mean["rebuilt"], rtol=0, atol=1e-12)
old_mean.select("compound_id", "endpoint", pl.col("prediction").alias("probability")).with_columns(
    pl.lit("historical_stack_mordred").alias("method"), pl.lit(cfg["split_id"]).alias("split_id")
).write_parquet(out / "test-probabilities.parquet")
auxiliary = {}
aux_memberships = []
heads = None
for source in cfg["auxiliary"]:
    summary = json.loads((inp / source["summary"]).read_text())
    excluded = tuple(source["excluded_folds"])
    assert summary["config"]["excluded_folds"] == list(excluded)
    membership = pl.read_parquet(inp / source["membership"]).join(
        assign.select("compound_id", "fold"), on="compound_id", validate="1:1"
    )
    expected = membership["fold"].is_in(excluded)
    np.testing.assert_array_equal(expected, membership["role"] == "outer_holdout")
    aux_memberships.append(
        dict(
            excluded_folds=list(excluded),
            run=source["run"],
            training=int((membership["role"] == "inner_train").sum()),
            validation=int((membership["role"] == "inner_validation").sum()),
        )
    )
    if heads is None:
        heads = summary["heads"]
    assert summary["heads"] == heads and len(heads) == 4
    predictions = pl.concat([pl.read_parquet(inp / name) for name in source["predictions"]])
    assert predictions["compound_id"].n_unique() == predictions.height == 6896
    auxiliary[excluded] = predictions.select("compound_id", *heads)
torch.set_num_threads(4)
rows, metrics, lineage = [], [], []
for outer in range(5):
    torch.manual_seed(42)
    inner = (outer + 1) % 5
    others = [fold for fold in range(5) if fold not in [outer, inner]]
    for endpoint in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        start = time.monotonic()
        frame = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(assign.select("compound_id", "fold"), on="compound_id", validate="1:1")
            .join(emb, on="compound_id", validate="1:1")
            .join(desc, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        ff = frame["fold"].to_numpy()
        train = ~np.isin(ff, [outer, inner])
        val = ff == inner
        y = frame["value"].to_numpy()
        raw = frame.select(ecols).to_numpy().astype(np.float32)
        pca = PCA(n_components=256, random_state=42 + outer)
        reduced = np.empty((frame.height, 256), dtype=np.float32)
        reduced[train] = pca.fit_transform(raw[train])
        reduced[val] = pca.transform(raw[val])
        assert pca.n_samples_ == train.sum()
        np.testing.assert_allclose(pca.mean_, raw[train].mean(0), atol=1e-6)
        aux = np.zeros((frame.height, 4), dtype=np.float32)
        for other in others:
            excluded = tuple(sorted([outer, inner, other]))
            values = (
                frame.select("compound_id")
                .join(auxiliary[excluded], on="compound_id", validate="1:1", maintain_order="left")
                .select(heads)
                .to_numpy()
            )
            aux[ff == other] = values[ff == other]
            aux[val] += values[val] / 3
            lineage.append(
                dict(
                    outer_fold=outer,
                    inner_fold=inner,
                    training_fold=other,
                    excluded_folds=list(excluded),
                    endpoint=endpoint,
                )
            )
        d = frame.select(dcols).to_numpy().astype(np.float64)
        d[~np.isfinite(d)] = np.nan
        keep = np.isfinite(d[train]).any(0)
        keep &= np.nan_to_num(np.nanstd(d[train], axis=0), nan=0) > 0
        imputer = SimpleImputer(strategy="median")
        dtrain = imputer.fit_transform(d[train][:, keep])
        dval = imputer.transform(d[val][:, keep])
        xtrain = np.column_stack([reduced[train], aux[train], dtrain])
        xval = np.column_stack([reduced[val], aux[val], dval])
        assert np.isfinite(xtrain).all() and np.isfinite(xval).all()
        name = f"{endpoint}-outer{outer}"
        joblib.dump(
            dict(
                pca=pca,
                imputer=imputer,
                descriptor_columns=[c for c, k in zip(dcols, keep, strict=True) if k],
                aux_heads=heads,
            ),
            out / f"{name}-preprocessors.joblib",
        )
        model = TabPFNClassifier.create_default_for_version(
            ModelVersion.V2_6,
            model_path=str(inp / cfg["weights"]),
            device="cuda",
            n_estimators=8,
            random_state=42 + outer,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
        )
        model.fit(xtrain, y[train])
        probability = model.predict_proba(xval)[:, list(model.classes_).index(1)]
        model.save_fit_state(out / f"{name}.tabpfn_fit")
        frame.select("compound_id", "fold").with_columns(
            pl.Series("training", train),
            pl.Series("inner_validation", val),
            pl.lit(outer).alias("outer_fold"),
        ).write_parquet(out / f"{name}-membership.parquet")
        rows.append(
            frame.filter(pl.Series(val))
            .select("compound_id", "endpoint", pl.col("value").alias("y_true"), "fold")
            .with_columns(
                pl.Series("probability", probability.astype(float)),
                pl.lit(outer).alias("outer_fold"),
                pl.lit(inner).alias("inner_fold"),
                pl.lit("historical_stack_mordred").alias("method"),
                pl.lit(cfg["split_id"]).alias("split_id"),
            )
        )
        metrics.append(
            dict(
                endpoint=endpoint,
                outer_fold=outer,
                inner_fold=inner,
                n_train=int(train.sum()),
                n_validation=int(val.sum()),
                dimensions=xtrain.shape[1],
                AUROC=roc_auc_score(y[val], probability),
                AP=average_precision_score(y[val], probability),
                seconds=time.monotonic() - start,
            )
        )
        print(json.dumps(metrics[-1]), flush=True)
        pl.concat(rows).write_parquet(out / "inner-probabilities.parquet")
        pl.DataFrame(metrics).write_csv(out / "inner-metrics.csv")
        del model
        gc.collect()
        torch.cuda.empty_cache()
assert pl.concat(rows).height == 5081
(out / "summary.json").write_text(
    json.dumps(
        dict(
            readout_fits=10,
            auxiliary_fits=10,
            rows=5081,
            feature_lineage=lineage,
            auxiliary_membership_audit=aux_memberships,
            outer_oof_unchanged=cfg["outer_oof_run"],
            test_unchanged=cfg["test_run"],
        ),
        indent=2,
    )
)
