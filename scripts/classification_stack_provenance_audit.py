"""Independent historical stack reconstruction lineage and preprocessing audit."""

import hashlib
import itertools
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.metrics import roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def p(i):
    return next(inp.glob(f"{i:03d}-*"))


def read(i):
    return pl.read_parquet(p(i))


labels = read(cfg["labels"])
assign = read(cfg["assignments"]).filter(pl.col("strategy") == "umap")
targets = read(cfg["targets"])
heads = sorted(targets["head"].unique())
wide = (
    targets.pivot(on="head", index="compound_id", values="value")
    .join(assign.select("compound_id", "fold"), on="compound_id", validate="1:1")
    .sort("compound_id")
)
y = wide.select(heads).to_numpy().astype(np.float32)
ff = wide["fold"].to_numpy()
auxrows = []
for source in cfg["auxiliary"]:
    summary = json.loads(p(source["summary"]).read_text())
    metadata = json.loads(p(source["metadata"]).read_text())
    c = summary["config"]
    excluded = source["excluded_folds"]
    assert c["excluded_folds"] == excluded and c["outer_fold"] == min(excluded)
    assert summary["heads"] == heads and c["seed"] == 42
    eligible = np.flatnonzero(~np.isin(ff, excluded))
    order = np.random.default_rng(42).permutation(eligible)
    nval = max(1, int(len(eligible) * 0.1))
    validation, training = order[:nval], order[nval:]
    expected = np.full(len(wide), "outer_holdout", dtype=object)
    expected[validation] = "inner_validation"
    expected[training] = "inner_train"
    membership = read(source["membership"]).sort("compound_id")
    assert membership["compound_id"].equals(wide["compound_id"])
    np.testing.assert_array_equal(membership["role"], expected)
    mean = np.nanmean(y[training], axis=0).astype(np.float32)
    std = np.nanstd(y[training], axis=0).astype(np.float32)
    std = np.where(std < 1e-6, 1, std)
    for doc in [summary, metadata]:
        np.testing.assert_allclose(doc["target_means"], mean, rtol=0, atol=1e-7)
        np.testing.assert_allclose(doc["target_stds"], std, rtol=0, atol=1e-7)
    history = [json.loads(line) for line in p(source["history"]).read_text().splitlines()]
    best = min(history, key=lambda r: r["val_mse_z"])
    checkpoint = torch.load(p(source["checkpoint"]), map_location="cpu", weights_only=False)
    assert best["epoch"] == summary["best_epoch"] == checkpoint["epoch"]
    np.testing.assert_allclose(
        best["val_mse_z"], checkpoint["validation_mse_z"], rtol=0, atol=1e-12
    )
    assert (
        hashlib.sha256(p(source["checkpoint"]).read_bytes()).hexdigest()
        == summary["best_checkpoint_sha256"]
    )
    predictions = pl.concat([read(i) for i in source["predictions"]])
    assert predictions.height == predictions["compound_id"].n_unique() == 6896
    assert np.isfinite(predictions.select(heads).to_numpy()).all()
    auxrows.append(
        dict(
            run=source["run"],
            excluded_folds=excluded,
            n_training=len(training),
            n_validation=len(validation),
            best_epoch=best["epoch"],
            epochs_run=len(history),
            scalers_exact=True,
            masks_exact=True,
        )
    )
assert {tuple(r["excluded_folds"]) for r in auxrows} == set(itertools.combinations(range(5), 3))
embedding = pl.concat([read(i) for i in cfg["embedding"]])
ecols = sorted(c for c in embedding.columns if c.startswith("chemeleon_"))
desc = pl.concat([read(i) for i in cfg["mordred"]])
dcols = [c for c in desc.columns if c != "compound_id"]
inner = read(cfg["inner"])
assert inner.height == 5081
metrics = pl.read_csv(p(cfg["metrics"]))
readout_checks = []
for item in cfg["readouts"]:
    ep, outer = item["endpoint"], item["outer_fold"]
    innerfold = (outer + 1) % 5
    frame = (
        labels.filter(pl.col("endpoint") == ep)
        .join(assign.select("compound_id", "fold"), on="compound_id", validate="1:1")
        .join(embedding, on="compound_id", validate="1:1")
        .join(desc, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    fold = frame["fold"].to_numpy()
    train = ~np.isin(fold, [outer, innerfold])
    val = fold == innerfold
    membership = read(item["membership"]).sort("compound_id")
    assert membership["compound_id"].equals(frame["compound_id"])
    np.testing.assert_array_equal(membership["training"], train)
    np.testing.assert_array_equal(membership["inner_validation"], val)
    pp = joblib.load(p(item["preprocessor"]))
    raw = frame.select(ecols).to_numpy().astype(np.float32)
    assert pp["pca"].n_samples_ == train.sum() and pp["pca"].n_components_ == 256
    np.testing.assert_allclose(pp["pca"].mean_, raw[train].mean(0), rtol=0, atol=1e-6)
    x = frame.select(dcols).to_numpy().astype(np.float64)
    x[~np.isfinite(x)] = np.nan
    keep = np.isfinite(x[train]).any(0) & (np.nan_to_num(np.nanstd(x[train], axis=0), nan=0) > 0)
    assert pp["descriptor_columns"] == [c for c, k in zip(dcols, keep, strict=True) if k]
    np.testing.assert_allclose(
        pp["imputer"].statistics_, np.nanmedian(x[train][:, keep], axis=0), rtol=0, atol=1e-12
    )
    assert pp["aux_heads"] == heads
    predicted = inner.filter((pl.col("endpoint") == ep) & (pl.col("outer_fold") == outer)).sort(
        "compound_id"
    )
    expected = frame.filter(pl.Series(val))
    assert predicted["compound_id"].equals(expected["compound_id"])
    np.testing.assert_array_equal(predicted["y_true"], expected["value"])
    assert (predicted["fold"] == innerfold).all() and predicted["probability"].is_between(
        0, 1
    ).all()
    mm = metrics.filter((pl.col("endpoint") == ep) & (pl.col("outer_fold") == outer))
    np.testing.assert_allclose(
        roc_auc_score(predicted["y_true"], predicted["probability"]),
        mm["AUROC"][0],
        atol=1e-12,
        rtol=0,
    )
    readout_checks.append(
        dict(
            endpoint=ep,
            outer_fold=outer,
            n_train=int(train.sum()),
            n_inner=int(val.sum()),
            PCA_train_only=True,
            descriptor_train_only=True,
        )
    )
summary = json.loads(p(cfg["summary"]).read_text())
expected_lineage = []
for outer in range(5):
    innerfold = (outer + 1) % 5
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        for fold in range(5):
            if fold not in [outer, innerfold]:
                expected_lineage.append(
                    (outer, innerfold, fold, ep, tuple(sorted([outer, innerfold, fold])))
                )
actual = [
    (
        r["outer_fold"],
        r["inner_fold"],
        r["training_fold"],
        r["endpoint"],
        tuple(r["excluded_folds"]),
    )
    for r in summary["feature_lineage"]
]
assert sorted(actual) == sorted(expected_lineage)
for original, new, folded in [
    (cfg["old_outer"], cfg["new_outer"], True),
    (cfg["old_test"], cfg["new_test"], True),
    (cfg["old_mean"], cfg["new_mean"], False),
]:
    keys = ["compound_id", "endpoint"] + (["fold"] if folded else [])
    a = (
        read(original)
        .filter(pl.col("endpoint").is_in(["CYP2D6_is_TDI", "CYP3A4_is_TDI"]))
        .sort(keys)
    )
    b = read(new).sort(keys)
    for key in keys:
        assert a[key].equals(b[key])
    np.testing.assert_array_equal(a["prediction"], b["probability"])
(out / "audit.json").write_text(
    json.dumps(
        dict(
            auxiliary=auxrows,
            readouts=readout_checks,
            inner_rows=5081,
            feature_lineage_rows=len(actual),
            old_outer_and_test_exactly_unchanged=True,
            scope="Independent saved artifacts and source review; no prediction refit",
        ),
        indent=2,
    )
)
print("PASS", len(auxrows), "auxiliary", len(readout_checks), "readouts", flush=True)
