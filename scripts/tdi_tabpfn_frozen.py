"""Frozen external representations for TDI classification on fixed UMAP folds."""

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
from sklearn.metrics import balanced_accuracy_score, log_loss, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion


def choose_threshold(y, p):
    candidates = np.linspace(0.01, 0.99, 99)
    return float(max(candidates, key=lambda t: (matthews_corrcoef(y, p >= t), -abs(t - 0.5), -t)))


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / cfg["labels"]).select("compound_id", "endpoint", "value")
    assigned = pl.read_parquet(inp / cfg["assignments"])
    assert assigned["split_id"].n_unique() == 1
    cfg["split_id"] = assigned["split_id"][0]
    mapping = pl.read_parquet(inp / cfg["test_mapping"])
    torch.set_num_threads(2)
    all_oof, all_test, metrics, thresholds = [], [], [], []
    endpoints = cfg.get("endpoints", ["CYP2D6_is_TDI", "CYP3A4_is_TDI"])
    for view in cfg["views"]:
        feature = pl.concat([pl.read_parquet(inp / name) for name in view["files"]])
        columns = [c for c in feature.columns if c.startswith(view["prefix"])]
        assert columns and feature["compound_id"].n_unique() == 6896
        for endpoint in endpoints:
            frame = (
                labels.filter(pl.col("endpoint") == endpoint)
                .join(assigned.select("compound_id", "fold"), on="compound_id", validate="1:1")
                .join(feature, on="compound_id", validate="1:1")
                .sort("compound_id")
            )
            assert frame.height == cfg["expected_counts"][endpoint]
            test = mapping.join(feature, on="compound_id", validate="1:1").sort("compound_id")
            assert test.height == 750
            x, xt = frame.select(columns).to_numpy(), test.select(columns).to_numpy()
            assert np.isfinite(x).all() and np.isfinite(xt).all()
            y = frame["value"].to_numpy().astype(int)
            assert set(y) == {0, 1} and set(frame["fold"]) == set(range(5))
            probs = np.full(len(y), np.nan)
            fold_test = []
            for fold in range(5):
                start = time.monotonic()
                val = frame["fold"].to_numpy() == fold
                scaler = StandardScaler().fit(x[~val])
                scaled_train = scaler.transform(x[~val])
                if view.get("preprocessing", "pca256") == "raw_scaled":
                    pca = None
                    train_x = scaled_train.astype(np.float32)
                    val_x = scaler.transform(x[val]).astype(np.float32)
                    test_x = scaler.transform(xt).astype(np.float32)
                else:
                    pca = PCA(n_components=256, random_state=42).fit(scaled_train)
                    train_x = pca.transform(scaled_train).astype(np.float32)
                    val_x = pca.transform(scaler.transform(x[val])).astype(np.float32)
                    test_x = pca.transform(scaler.transform(xt)).astype(np.float32)
                name = f"{view['name']}-{endpoint}-fold{fold}"
                joblib.dump(
                    dict(scaler=scaler, pca=pca, columns=columns), out / f"{name}-pca.joblib"
                )
                model = TabPFNClassifier.create_default_for_version(
                    ModelVersion.V2_6,
                    n_estimators=4,
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
                model.fit(train_x, y[~val])
                positive = list(model.classes_).index(1)
                probs[val] = model.predict_proba(val_x)[:, positive]
                tp = model.predict_proba(test_x)[:, positive].astype(float)
                model.save_fit_state(out / f"{name}.tabpfn_fit")
                frame.select("compound_id", "fold").with_columns(
                    pl.Series("training", ~val), pl.Series("validation", val)
                ).write_parquet(out / f"{name}-membership.parquet")
                fold_test.append(
                    test.select("compound_id").with_columns(
                        pl.lit(endpoint).alias("endpoint"),
                        pl.lit(view["name"]).alias("method"),
                        pl.lit(fold).alias("fold"),
                        pl.Series("probability", tp),
                        pl.lit(cfg["split_id"]).alias("split_id"),
                    )
                )
                print(
                    json.dumps(
                        dict(
                            method=view["name"],
                            endpoint=endpoint,
                            fold=fold,
                            train=int((~val).sum()),
                            validation=int(val.sum()),
                            seconds=time.monotonic() - start,
                        )
                    ),
                    flush=True,
                )
                del model
                gc.collect()
                torch.cuda.empty_cache()
            assert np.isfinite(probs).all()
            binary = np.zeros(len(y), dtype=int)
            chosen = np.zeros(len(y))
            for fold in range(5):
                val = frame["fold"].to_numpy() == fold
                threshold = choose_threshold(y[~val], probs[~val])
                chosen[val] = threshold
                binary[val] = probs[val] >= threshold
                thresholds.append(
                    dict(
                        method=view["name"],
                        endpoint=endpoint,
                        fold=fold,
                        threshold=threshold,
                        purpose="other_fold_OOF_MCC",
                    )
                )
            deployment = choose_threshold(y, probs)
            thresholds.append(
                dict(
                    method=view["name"],
                    endpoint=endpoint,
                    fold=-1,
                    threshold=deployment,
                    purpose="all_OOF_test_only",
                )
            )
            oof = frame.select("compound_id", "endpoint", "fold").with_columns(
                pl.Series("y_true", y),
                pl.Series("probability", probs),
                pl.Series("prediction_05", (probs >= 0.5).astype(int)),
                pl.Series("prediction_crossfit", binary),
                pl.Series("threshold", chosen),
                pl.lit(view["name"]).alias("method"),
                pl.lit(cfg["split_id"]).alias("split_id"),
            )
            all_oof.append(oof)
            all_test.extend(fold_test)
            metrics.append(
                dict(
                    method=view["name"],
                    endpoint=endpoint,
                    n=len(y),
                    positives=int(y.sum()),
                    mcc_05=matthews_corrcoef(y, probs >= 0.5),
                    mcc_crossfit=matthews_corrcoef(y, binary),
                    roc_auc=roc_auc_score(y, probs),
                    log_loss=log_loss(y, probs),
                    balanced_accuracy_05=balanced_accuracy_score(y, probs >= 0.5),
                    deployment_threshold=deployment,
                )
            )
            pl.concat(all_oof).write_parquet(out / "oof-probabilities.parquet")
            pl.concat(all_test).write_parquet(out / "test-fold-probabilities.parquet")
            pl.DataFrame(metrics).write_csv(out / "metrics.csv")
            pl.DataFrame(thresholds).write_csv(out / "thresholds.csv")
    mean = (
        pl.concat(all_test)
        .group_by("compound_id", "endpoint", "method", "split_id")
        .agg(pl.col("probability").mean(), pl.len().alias("n_folds"))
    )
    expected_test_rows = 750 * len(cfg["views"]) * len(endpoints)
    assert (mean["n_folds"] == 5).all() and mean.height == expected_test_rows
    deployment = (
        pl.DataFrame(thresholds)
        .filter(pl.col("fold") == -1)
        .select("method", "endpoint", "threshold")
    )
    mean = mean.join(deployment, on=["method", "endpoint"], validate="m:1").with_columns(
        (pl.col("probability") >= 0.5).cast(pl.Int8).alias("prediction_05"),
        (pl.col("probability") >= pl.col("threshold"))
        .cast(pl.Int8)
        .alias("prediction_oof_threshold"),
    )
    mean.write_parquet(out / "test-probabilities.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                fits=5 * len(cfg["views"]) * len(endpoints),
                config=cfg,
                threshold_note="Other-fold OOF thresholds; overlapping fits, not fully nested.",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
