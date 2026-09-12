"""Independent inner UMAP threshold fits applied to frozen outer probabilities."""

import gc
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.decomposition import PCA
from sklearn.metrics import matthews_corrcoef
from sklearn.preprocessing import StandardScaler
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion
from tdi_tabpfn_frozen import choose_threshold

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
labels = pl.read_parquet(inp / cfg["labels"]).select("compound_id", "endpoint", "value")
assigned = pl.read_parquet(inp / cfg["assignments"])
outer = pl.read_parquet(inp / cfg["outer_oof"])
test = pl.read_parquet(inp / cfg["outer_test"])
torch.set_num_threads(2)
results, thresholds, validations = [], [], []
for view in cfg["views"]:
    feature = pl.concat([pl.read_parquet(inp / name) for name in view["files"]])
    columns = [c for c in feature.columns if c.startswith(view["prefix"])]
    for endpoint in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        frame = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(assigned.select("compound_id", "fold"), on="compound_id", validate="1:1")
            .join(feature, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x, y, folds = (
            frame.select(columns).to_numpy(),
            frame["value"].to_numpy().astype(int),
            frame["fold"].to_numpy(),
        )
        for fold in range(5):
            inner = (fold + 1) % 5
            train, val = (folds != fold) & (folds != inner), folds == inner
            assert not (train & (folds == fold)).any()
            scaler = StandardScaler().fit(x[train])
            pca = PCA(n_components=256, random_state=42).fit(scaler.transform(x[train]))
            name = f"{view['name']}-{endpoint}-outer{fold}"
            joblib.dump(dict(scaler=scaler, pca=pca, columns=columns), out / f"{name}-pca.joblib")
            model = TabPFNClassifier.create_default_for_version(
                ModelVersion.V2_6,
                n_estimators=4,
                model_path=str(inp / cfg["weights"]),
                device="cuda",
                random_state=42 + fold,
                n_preprocessing_jobs=1,
                memory_saving_mode="auto",
                show_progress_bar=False,
                ignore_pretraining_limits=True,
                fit_mode="fit_preprocessors",
            )
            model.fit(pca.transform(scaler.transform(x[train])).astype(np.float32), y[train])
            probs = model.predict_proba(pca.transform(scaler.transform(x[val])).astype(np.float32))[
                :, list(model.classes_).index(1)
            ]
            model.save_fit_state(out / f"{name}.tabpfn_fit")
            threshold = choose_threshold(y[val], probs)
            thresholds.append(
                dict(
                    method=view["name"],
                    endpoint=endpoint,
                    fold=fold,
                    inner_fold=inner,
                    threshold=threshold,
                )
            )
            current = outer.filter(
                (pl.col("method") == view["name"])
                & (pl.col("endpoint") == endpoint)
                & (pl.col("fold") == fold)
            )
            assert set(current["compound_id"]) == set(
                frame.filter(pl.col("fold") == fold)["compound_id"]
            )
            results.append(
                current.drop("threshold").with_columns(
                    pl.lit(threshold).alias("inner_threshold"),
                    (pl.col("probability") >= threshold).cast(pl.Int8).alias("prediction_inner"),
                )
            )
            validations.append(
                frame.filter(pl.Series(val))
                .select("compound_id", "endpoint")
                .with_columns(
                    pl.Series("y_true", y[val]),
                    pl.Series("probability", probs),
                    pl.lit(fold).alias("outer_fold"),
                    pl.lit(inner).alias("inner_fold"),
                    pl.lit(view["name"]).alias("method"),
                )
            )
            frame.select("compound_id", "fold").with_columns(
                pl.Series("training", train),
                pl.Series("inner_validation", val),
                pl.lit(fold).alias("outer_fold"),
            ).write_parquet(out / f"{name}-membership.parquet")
            print(json.dumps(thresholds[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
            pl.concat(results).write_parquet(out / "oof-inner-threshold.parquet")
            pl.DataFrame(thresholds).write_csv(out / "thresholds.csv")
pl.concat(validations).write_parquet(out / "inner-validation-probabilities.parquet")
rows = []
for (method, endpoint), table in pl.concat(results).group_by("method", "endpoint"):
    rows.append(
        dict(
            method=method,
            endpoint=endpoint,
            n=table.height,
            mcc_05=matthews_corrcoef(table["y_true"], table["prediction_05"]),
            mcc_inner=matthews_corrcoef(table["y_true"], table["prediction_inner"]),
        )
    )
pl.DataFrame(rows).write_csv(out / "metrics.csv")
test.join(pl.DataFrame(thresholds), on=["method", "endpoint", "fold"], validate="m:1").with_columns(
    (pl.col("probability") >= pl.col("threshold")).cast(pl.Int8).alias("prediction_inner")
).write_parquet(out / "test-fold-inner-threshold.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            fits=30,
            threshold_policy="Train excludes outer and next fold; next fold selects threshold",
            source_run=cfg["outer_run"],
        ),
        indent=2,
    )
)
