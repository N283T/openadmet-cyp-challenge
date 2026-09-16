"""Fixed-readout screen for generic KERMT and MoLFormer CYP2D6 signal."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from catboost import CatBoostRegressor
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from threadpoolctl import threadpool_limits
from xgboost import XGBRegressor

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]


def read(index: int) -> pl.DataFrame:
    paths = list(inp.glob(f"{index:03d}-*"))
    if len(paths) != 1:
        raise ValueError(f"Expected one input at index {index}, found {paths}")
    return pl.read_parquet(paths[0])


def st_rae(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    y = frame["y_true"].to_numpy()
    lower = frame["conf_low"].to_numpy()
    upper = frame["conf_high"].to_numpy()
    error = np.maximum(lower - prediction, 0) + np.maximum(prediction - upper, 0)
    center = np.full_like(y, y.mean())
    denominator = (np.maximum(lower - center, 0) + np.maximum(center - upper, 0)).sum()
    return float(error.sum() / denominator)


def make_model(method: str, dimensions: int, seed: int) -> object:
    if method == "pca64_ridge":
        return make_pipeline(
            StandardScaler(),
            PCA(n_components=min(64, dimensions), svd_solver="randomized", random_state=seed),
            StandardScaler(),
            Ridge(alpha=10.0),
        )
    if method == "pca64_rbf_svr":
        return make_pipeline(
            StandardScaler(),
            PCA(n_components=min(64, dimensions), svd_solver="randomized", random_state=seed),
            StandardScaler(),
            SVR(C=10.0, epsilon=0.1, gamma="scale"),
        )
    if method == "pls16":
        return make_pipeline(
            StandardScaler(), PLSRegression(n_components=min(16, dimensions), scale=False)
        )
    if method == "extra_trees":
        return ExtraTreesRegressor(
            n_estimators=cfg["tree_estimators"],
            min_samples_leaf=3,
            max_features=0.5,
            n_jobs=cfg["threads"],
            random_state=seed,
        )
    if method.startswith("extra_trees_l"):
        leaf = int(method.rsplit("l", maxsplit=1)[1])
        return ExtraTreesRegressor(
            n_estimators=cfg["tree_estimators"],
            min_samples_leaf=leaf,
            max_features=0.5,
            n_jobs=cfg["threads"],
            random_state=seed,
        )
    if method == "extra_trees_f025":
        return ExtraTreesRegressor(
            n_estimators=cfg["tree_estimators"],
            min_samples_leaf=3,
            max_features=0.25,
            n_jobs=cfg["threads"],
            random_state=seed,
        )
    if method == "random_forest":
        return RandomForestRegressor(
            n_estimators=cfg["tree_estimators"],
            min_samples_leaf=3,
            max_features=0.5,
            n_jobs=cfg["threads"],
            random_state=seed,
        )
    if method == "catboost_mae":
        return CatBoostRegressor(
            iterations=500,
            depth=6,
            learning_rate=0.03,
            l2_leaf_reg=5.0,
            loss_function="MAE",
            random_seed=seed,
            thread_count=cfg["threads"],
            verbose=False,
            allow_writing_files=False,
        )
    if method == "xgboost":
        return XGBRegressor(
            n_estimators=500,
            max_depth=3,
            learning_rate=0.03,
            min_child_weight=5,
            subsample=0.8,
            colsample_bytree=0.5,
            reg_lambda=5.0,
            objective="reg:squarederror",
            tree_method="hist",
            n_jobs=cfg["threads"],
            random_state=seed,
        )
    raise ValueError(method)


if cfg.get("feature_mode") == "molformer_pooling":
    assignments, labels, scores, nested, latest, current, pooled = [
        read(index) for index in range(7)
    ]
    score_columns = sorted(column for column in scores.columns if column.startswith("score_"))
    pooled_columns = sorted(
        column for column in pooled.columns if column.startswith("molformer_cls_mean_")
    )
    if (len(score_columns), len(pooled_columns)) != (47, 1536):
        raise ValueError("Unexpected pooling feature dimensions")
    cls_columns = pooled_columns[:768]
    mean_columns = pooled_columns[768:]
    static = scores.select("compound_id", *score_columns).join(
        pooled, on="compound_id", validate="1:1"
    )
    views = {
        "molformer_mean": mean_columns,
        "molformer_cls_mean": cls_columns + mean_columns,
        "scores_molformer_mean": score_columns + mean_columns,
        "scores_molformer_cls_mean": score_columns + cls_columns + mean_columns,
    }
    if cfg.get("feature_views"):
        views = {name: views[name] for name in cfg["feature_views"]}
else:
    assignments, labels, scores, nested, latest, current, kermt, molformer = [
        read(index) for index in range(8)
    ]
    score_columns = sorted(column for column in scores.columns if column.startswith("score_"))
    kermt_columns = sorted(column for column in kermt.columns if column.startswith("kermt_"))
    molformer_columns = sorted(
        column for column in molformer.columns if column.startswith("molformer_")
    )
    if (len(score_columns), len(kermt_columns), len(molformer_columns)) != (47, 3200, 768):
        raise ValueError("Unexpected feature dimensions")
    static = (
        scores.select("compound_id", *score_columns)
        .join(kermt, on="compound_id", validate="1:1")
        .join(molformer, on="compound_id", validate="1:1")
    )
    views = {
        "kermt": kermt_columns,
        "molformer": molformer_columns,
        "scores_kermt": score_columns + kermt_columns,
        "scores_molformer": score_columns + molformer_columns,
    }
labels = labels.filter(pl.col("endpoint") == cfg["endpoint"])
base = labels.join(assignments, on="compound_id", validate="1:1").sort("compound_id")
latest = latest.filter(pl.col("endpoint") == cfg["endpoint"]).sort(keys)
current = current.filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("component_count") == 2)
).sort(keys)
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not latest[column].equals(current[column]):
        raise ValueError(f"Current/latest mismatch in {column}")

rows = []
fit_rows = []
memberships = []
threadpool_limits(limits=cfg["threads"])
for fold in range(5):
    data = (
        base.join(static, on="compound_id", validate="1:1")
        .join(
            nested.filter(pl.col("outer_fold") == fold).select(
                "compound_id", *[f"aux_{index}" for index in range(4)]
            ),
            on="compound_id",
            validate="1:1",
        )
        .sort("compound_id")
    )
    validation = (data["fold"].to_numpy() == fold) & ~data["is_anchor"].to_numpy()
    train = ~validation
    y = data["value"].to_numpy().astype(float)
    memberships.append(
        data.select("compound_id", "endpoint", "fold", "group", "is_anchor").with_columns(
            pl.lit(fold).alias("outer_fold"), pl.Series("validation", validation)
        )
    )
    for view, base_columns in views.items():
        columns = base_columns + [f"aux_{index}" for index in range(4)]
        x = data.select(columns).to_numpy().astype(np.float32)
        for method in cfg["methods"]:
            model = clone(make_model(method, x.shape[1], cfg["seed"] + fold))
            model.fit(x[train], y[train])
            prediction = np.asarray(model.predict(x[validation])).reshape(-1).astype(float)
            directory = out / f"fold{fold}-{view}-{method}"
            directory.mkdir()
            joblib.dump({"model": model, "columns": columns}, directory / "model.joblib")
            held = data.filter(pl.Series(validation))
            row = held.select(
                "compound_id",
                "endpoint",
                "fold",
                "group",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            ).with_columns(
                pl.lit(view).alias("view"),
                pl.lit(method).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.Series("prediction", prediction),
            )
            rows.append(row)
            record = {
                "fold": fold,
                "view": view,
                "method": method,
                "train": int(train.sum()),
                "validation": int(validation.sum()),
                "dimensions": x.shape[1],
                "fold_ST_RAE": st_rae(row, prediction),
            }
            fit_rows.append(record)
            (directory / "fit.json").write_text(json.dumps(record, indent=2) + "\n")
            print(json.dumps(record), flush=True)
            pl.concat(rows).write_parquet(out / "oof-progress.parquet")

oof = pl.concat(rows).sort("view", "method", "fold", "compound_id")
expected = latest.height * len(views) * len(cfg["methods"])
if oof.height != expected:
    raise ValueError(f"Expected {expected} OOF rows, found {oof.height}")
oof.write_parquet(out / "oof.parquet")
pl.concat(memberships).write_parquet(out / "memberships.parquet")

baseline = latest.select(*keys, "prediction").rename({"prediction": "latest_prediction"})
metric_rows = []
diversity_rows = []
for (view, method), frame in oof.group_by("view", "method"):
    frame = frame.sort(keys).join(baseline, on=keys, validate="1:1")
    prediction = frame["prediction"].to_numpy()
    latest_prediction = frame["latest_prediction"].to_numpy()
    y = frame["y_true"].to_numpy()
    matched = (prediction - prediction.mean()) / prediction.std() * latest_prediction.std()
    matched += latest_prediction.mean()
    diversity_rows.append(
        {
            "view": view,
            "method": method,
            "prediction_correlation_latest": float(
                np.corrcoef(prediction, latest_prediction)[0, 1]
            ),
            "residual_correlation_latest": float(
                np.corrcoef(y - prediction, y - latest_prediction)[0, 1]
            ),
        }
    )
    for placement, candidate in (("raw", prediction), ("match_latest_moments", matched)):
        for weight in cfg["blend_weights"]:
            blended = (1 - weight) * latest_prediction + weight * candidate
            metric_rows.append(
                {
                    "view": view,
                    "method": method,
                    "placement": placement,
                    "new_weight": weight,
                    "ST_RAE": st_rae(frame, blended),
                    "MAE": float(np.abs(y - blended).mean()),
                    "Spearman": float(spearmanr(y, blended).statistic),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("ST_RAE", "view", "method")
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(diversity_rows).sort("residual_correlation_latest").write_csv(out / "diversity.csv")
pl.DataFrame(fit_rows).sort("fold", "view", "method").write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "feature_dimensions": {view: len(columns) + 4 for view, columns in views.items()},
            "baseline_latest_ST_RAE": st_rae(latest, latest["prediction"].to_numpy()),
            "baseline_current_ST_RAE": st_rae(current, current["prediction"].to_numpy()),
            "selection_note": (
                "Fixed readouts and blend weights; selected rows need bootstrap audit."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(metrics.head(30), flush=True)
