"""Cross-fold stacks and fixed combinations of retained CYP2D6 OOF signals."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

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


def select(frame: pl.DataFrame, view: str, method: str, name: str) -> pl.DataFrame:
    return (
        frame.filter((pl.col("view") == view) & (pl.col("method") == method))
        .sort(keys)
        .select(*keys, pl.col("prediction").alias(name))
    )


latest, current, molformer_oof, family_oof, ordinal_oof, pairwise_oof = [
    read(index) for index in range(6)
]
latest = latest.filter(pl.col("endpoint") == cfg["endpoint"]).sort(keys)
current = current.filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("component_count") == 2)
).sort(keys)
molformer = molformer_oof.filter(
    (pl.col("view") == "scores_molformer_cls_mean")
    & (pl.col("method") == "random_forest")
).sort(keys)
base = molformer.select(
    *keys, "group", "is_neighbor", "y_true", "conf_low", "conf_high"
).join(
    latest.select(*keys, pl.col("prediction").alias("latest")), on=keys, validate="1:1"
)
components = [
    select(molformer_oof, "scores_molformer_cls_mean", "random_forest", "molformer_rf"),
    select(family_oof, "hybrid2699", "catboost", "hybrid_catboost"),
    select(ordinal_oof, "scores51", "ordinal_logit5", "ordinal5"),
    select(pairwise_oof, "scores51", "xgb_pairwise", "xgb_pairwise"),
]
for component in components:
    base = base.join(component, on=keys, validate="1:1")
if base.height != cfg["expected_rows"]:
    raise ValueError("OOF coverage mismatch")
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not latest[column].equals(current[column]) or not latest[column].equals(base[column]):
        raise ValueError(f"Baseline mismatch in {column}")

y = base["y_true"].to_numpy()
latest_prediction = base["latest"].to_numpy()
raw_columns = ["latest", "molformer_rf", "hybrid_catboost", "ordinal5", "xgb_pairwise"]
x = base.select(raw_columns).to_numpy()

matched_molformer = (
    (base["molformer_rf"].to_numpy() - base["molformer_rf"].mean())
    / base["molformer_rf"].std()
    * base["latest"].std()
    + base["latest"].mean()
)
matched_catboost = (
    (base["hybrid_catboost"].to_numpy() - base["hybrid_catboost"].mean())
    / base["hybrid_catboost"].std()
    * base["latest"].std()
    + base["latest"].mean()
)
retained = {
    "molformer_rf25": 0.75 * latest_prediction + 0.25 * matched_molformer,
    "hybrid_catboost50": 0.5 * latest_prediction + 0.5 * matched_catboost,
    "ordinal5_10": 0.9 * latest_prediction + 0.1 * base["ordinal5"].to_numpy(),
    "xgb_pairwise10": 0.9 * latest_prediction + 0.1 * base["xgb_pairwise"].to_numpy(),
}
predictions = {
    "fixed_molformer_catboost": 0.5
    * (retained["molformer_rf25"] + retained["hybrid_catboost50"]),
    "fixed_molformer_ordinal": 0.5
    * (retained["molformer_rf25"] + retained["ordinal5_10"]),
    "fixed_molformer_pairwise": 0.5
    * (retained["molformer_rf25"] + retained["xgb_pairwise10"]),
    "fixed_all_four": np.mean(np.vstack(list(retained.values())), axis=0),
}

fold_ids = base["fold"].to_numpy()
stack_methods = {
    "ridge_stack_a1": lambda: make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
    "ridge_stack_a10": lambda: make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
    "ridge_stack_a100": lambda: make_pipeline(StandardScaler(), Ridge(alpha=100.0)),
    "positive_linear_stack": lambda: LinearRegression(positive=True),
}
fit_rows = []
for method, factory in stack_methods.items():
    prediction = np.empty(base.height)
    for fold in range(5):
        validation = fold_ids == fold
        train = ~validation
        model = factory()
        model.fit(x[train], y[train])
        prediction[validation] = model.predict(x[validation])
        directory = out / f"fold{fold}-{method}"
        directory.mkdir()
        joblib.dump({"model": model, "columns": raw_columns}, directory / "model.joblib")
        fit_rows.append(
            {
                "method": method,
                "fold": fold,
                "train": int(train.sum()),
                "validation": int(validation.sum()),
            }
        )
    predictions[method] = prediction

oof_rows = []
metric_rows = []
fold_rows = []
bootstrap_rows = []
rng = np.random.default_rng(cfg["bootstrap_seed"])
group_to_indices = {
    group: indices["index"].to_numpy()
    for (group,), indices in base.with_row_index("index").group_by("group")
}
groups = np.asarray(sorted(group_to_indices))
baseline_score = st_rae(base, latest_prediction)
baseline_spearman = float(spearmanr(y, latest_prediction).statistic)
for method, prediction in predictions.items():
    oof_rows.append(
        base.select(*keys, "group", "is_neighbor", "y_true", "conf_low", "conf_high")
        .with_columns(pl.lit(method).alias("method"), pl.Series("prediction", prediction))
    )
    score = st_rae(base, prediction)
    metric_rows.append(
        {
            "method": method,
            "ST_RAE": score,
            "delta_ST_RAE": score - baseline_score,
            "Spearman": float(spearmanr(y, prediction).statistic),
            "delta_Spearman": float(spearmanr(y, prediction).statistic) - baseline_spearman,
            "prediction_correlation_latest": float(
                np.corrcoef(prediction, latest_prediction)[0, 1]
            ),
        }
    )
    for fold in range(5):
        mask = fold_ids == fold
        fold_frame = base.filter(pl.Series(mask))
        fold_rows.append(
            {
                "method": method,
                "fold": fold,
                "baseline_ST_RAE": st_rae(fold_frame, latest_prediction[mask]),
                "candidate_ST_RAE": st_rae(fold_frame, prediction[mask]),
            }
        )
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled_groups = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([group_to_indices[group] for group in sampled_groups])
        sampled = base[indices]
        deltas.append(
            st_rae(sampled, prediction[indices])
            - st_rae(sampled, latest_prediction[indices])
        )
    values = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "method": method,
            "delta_q025": float(np.quantile(values, 0.025)),
            "delta_q50": float(np.quantile(values, 0.5)),
            "delta_q975": float(np.quantile(values, 0.975)),
            "probability_improvement": float(np.mean(values < 0)),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("ST_RAE")
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(fold_rows).with_columns(
    (pl.col("candidate_ST_RAE") - pl.col("baseline_ST_RAE")).alias("delta_ST_RAE")
).sort("method", "fold").write_csv(out / "folds.csv")
pl.DataFrame(bootstrap_rows).sort("delta_q50").write_csv(out / "bootstrap.csv")
pl.DataFrame(fit_rows).sort("method", "fold").write_csv(out / "fits.csv")
pl.concat(oof_rows).sort("method", *keys).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "baseline_ST_RAE": baseline_score,
            "baseline_Spearman": baseline_spearman,
            "note": "Meta-stacks fit on four OOF folds and predict the excluded fifth fold.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
print(pl.DataFrame(bootstrap_rows).sort("delta_q50"), flush=True)
