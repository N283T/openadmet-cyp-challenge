"""Jeremy600 frozen features with matched TabICL8 and TabPFN8 readouts."""

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
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler
from tabicl import TabICLClassifier
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion


def choose_threshold(y, p):
    return float(
        max(
            np.linspace(0.01, 0.99, 99),
            key=lambda t: (matthews_corrcoef(y, p >= t), -abs(t - 0.5), -t),
        )
    )


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / cfg["labels"]).select("compound_id", "endpoint", "value")
    assigned = pl.read_parquet(inp / cfg["assignments"])
    mapping = pl.read_parquet(inp / cfg["test_mapping"])
    assert assigned["split_id"].n_unique() == 1
    split = assigned["split_id"][0]
    torch.set_num_threads(2)
    oo, ii, tt, metrics, thresholds = [], [], [], [], []
    for view in cfg["views"]:
        features = pl.concat([pl.read_parquet(inp / f) for f in view["files"]])
        columns = [c for c in features.columns if c.startswith(view["prefix"])]
        assert features.height == features["compound_id"].n_unique() == 6896
        for endpoint in cfg["expected_counts"]:
            frame = (
                labels.filter(pl.col("endpoint") == endpoint)
                .join(assigned.select("compound_id", "fold"), on="compound_id", validate="1:1")
                .join(features, on="compound_id", validate="1:1")
                .sort("compound_id")
            )
            test = mapping.join(features, on="compound_id", validate="1:1").sort("compound_id")
            assert frame.height == cfg["expected_counts"][endpoint] and test.height == 750
            x, xt = frame.select(columns).to_numpy(), test.select(columns).to_numpy()
            y, folds = frame["value"].to_numpy().astype(int), frame["fold"].to_numpy()
            assert np.isfinite(x).all() and np.isfinite(xt).all()
            endpoint_oof = []
            for fold in range(5):
                inner = (fold + 1) % 5
                pair = {}
                for phase in ["outer", "inner"]:
                    started = time.monotonic()
                    train = folds != fold
                    if phase == "inner":
                        train &= folds != inner
                    val = folds == (fold if phase == "outer" else inner)
                    assert not (train & val).any() and not (train & (folds == fold)).any()
                    scaler = StandardScaler().fit(x[train])
                    pca = PCA(n_components=256, random_state=42).fit(scaler.transform(x[train]))
                    tx = pca.transform(scaler.transform(x[train])).astype(np.float32)
                    vx = pca.transform(scaler.transform(x[val])).astype(np.float32)
                    name = f"{view['name']}-{endpoint}-fold{fold}-{phase}"
                    joblib.dump(
                        dict(scaler=scaler, pca=pca, columns=columns), out / f"{name}-pca.joblib"
                    )
                    frame.select("compound_id", "fold").with_columns(
                        pl.Series("training", train),
                        pl.Series("validation", val),
                        pl.lit(fold).alias("outer_fold"),
                        pl.lit(inner).alias("inner_fold"),
                    ).write_parquet(out / f"{name}-membership.parquet")
                    if view["backend"] == "tabicl":
                        model = TabICLClassifier(
                            n_estimators=8,
                            model_path=str(inp / cfg["weights"]),
                            allow_auto_download=False,
                            device="cuda",
                            random_state=42 + fold,
                            n_jobs=2,
                            batch_size=8,
                        )
                    else:
                        model = TabPFNClassifier.create_default_for_version(
                            ModelVersion.V2_6,
                            n_estimators=8,
                            model_path=str(inp / cfg["pfn_weights"]),
                            device="cuda",
                            random_state=42 + fold,
                            auto_scale_n_estimators=False,
                            n_preprocessing_jobs=1,
                            memory_saving_mode="auto",
                            show_progress_bar=False,
                            ignore_pretraining_limits=True,
                            fit_mode="fit_preprocessors",
                        )
                    model.fit(tx, y[train])
                    positive = list(model.classes_).index(1)
                    pred = model.predict_proba(vx)[:, positive].astype(float)
                    assert np.isfinite(pred).all() and ((pred >= 0) & (pred <= 1)).all()
                    table = (
                        frame.filter(pl.Series(val))
                        .select("compound_id", "endpoint", "fold")
                        .with_columns(
                            pl.Series("y_true", y[val]),
                            pl.Series("probability", pred),
                            pl.lit(view["name"]).alias("method"),
                            pl.lit(split).alias("split_id"),
                        )
                    )
                    pair[phase] = table
                    if phase == "outer":
                        test_x = pca.transform(scaler.transform(xt)).astype(np.float32)
                        tp = model.predict_proba(test_x)[:, positive].astype(float)
                        tt.append(
                            test.select("compound_id").with_columns(
                                pl.lit(endpoint).alias("endpoint"),
                                pl.lit(view["name"]).alias("method"),
                                pl.lit(fold).alias("fold"),
                                pl.lit(split).alias("split_id"),
                                pl.Series("probability", tp),
                            )
                        )
                    else:
                        ii.append(
                            table.with_columns(
                                pl.lit(fold).alias("outer_fold"), pl.lit(inner).alias("inner_fold")
                            )
                        )
                    if view["backend"] == "tabicl":
                        model.save(
                            out / f"{name}.tabicl_fit",
                            save_model_weights=False,
                            save_training_data=True,
                            save_kv_cache=False,
                        )
                    else:
                        model.save_fit_state(out / f"{name}.tabpfn_fit")
                    print(
                        json.dumps(
                            dict(
                                method=view["name"],
                                endpoint=endpoint,
                                fold=fold,
                                phase=phase,
                                train=int(train.sum()),
                                validation=int(val.sum()),
                                seconds=time.monotonic() - started,
                            )
                        ),
                        flush=True,
                    )
                    del model
                    gc.collect()
                    torch.cuda.empty_cache()
                threshold = choose_threshold(pair["inner"]["y_true"], pair["inner"]["probability"])
                thresholds.append(
                    dict(
                        method=view["name"],
                        endpoint=endpoint,
                        fold=fold,
                        inner_fold=inner,
                        threshold=threshold,
                    )
                )
                endpoint_oof.append(
                    pair["outer"].with_columns(
                        pl.lit(threshold).alias("inner_threshold"),
                        (pl.col("probability") >= 0.5).cast(pl.Int8).alias("prediction_05"),
                        (pl.col("probability") >= threshold)
                        .cast(pl.Int8)
                        .alias("prediction_inner"),
                    )
                )
            table = pl.concat(endpoint_oof)
            oo.append(table)
            metrics.append(
                dict(
                    method=view["name"],
                    endpoint=endpoint,
                    n=table.height,
                    roc_auc=roc_auc_score(table["y_true"], table["probability"]),
                    average_precision=average_precision_score(
                        table["y_true"], table["probability"]
                    ),
                    mcc_05=matthews_corrcoef(table["y_true"], table["prediction_05"]),
                    mcc_inner=matthews_corrcoef(table["y_true"], table["prediction_inner"]),
                )
            )
            pl.concat(oo).write_parquet(out / "outer-oof.parquet")
            pl.concat(ii).write_parquet(out / "inner-probabilities.parquet")
            pl.concat(tt).write_parquet(out / "fold-test-probabilities.parquet")
            pl.DataFrame(metrics).write_csv(out / "metrics.csv")
            pl.DataFrame(thresholds).write_csv(out / "thresholds.csv")
    assert pl.concat(oo).height == pl.concat(ii).height == 10162
    assert pl.concat(tt).height == 15000
    mean = (
        pl.concat(tt)
        .group_by("compound_id", "endpoint", "method", "split_id")
        .agg(pl.col("probability").mean(), pl.len().alias("n_folds"))
    )
    assert mean.height == 3000 and (mean["n_folds"] == 5).all()
    mean.write_parquet(out / "test-probabilities.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                fits=40,
                config=cfg,
                threshold_policy="Independent next-fold inner labels; same 0.01..0.99 grid as PFN",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
