"""Nested anchor-preserving stacking on a saved OOF prediction bank."""

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import polars as pl
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def fixed_calls(values, folds, rate):
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y, calls, prevalence):
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    tpr, fpr = float(calls[positive].mean()), float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def score(y, prediction, folds, settings):
    return expected_mcc(
        y,
        fixed_calls(prediction, folds, settings["call_rate"]),
        settings["target_prior"],
    )


def fit_anchor_simplex(x, y, regularization):
    count = x.shape[1]
    positive, negative = x[y == 1], x[y == 0]
    differences = (positive[:, None, :] - negative[None, :, :]).reshape(-1, count)
    temperature = cfg["pairwise_temperature"]
    anchor = np.zeros(count)
    anchor[0] = 1

    def objective(weights):
        margins = differences @ weights / temperature
        loss = float(np.logaddexp(0, -margins).mean())
        return loss + regularization * float(np.sum((weights - anchor) ** 2))

    def gradient(weights):
        margins = differences @ weights / temperature
        rank_gradient = (-expit(-margins)[:, None] * differences / temperature).mean(axis=0)
        return rank_gradient + 2 * regularization * (weights - anchor)

    result = minimize(
        objective,
        anchor,
        jac=gradient,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * count,
        constraints={"type": "eq", "fun": lambda weights: weights.sum() - 1.0},
        options={"maxiter": 500, "ftol": 1e-10},
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = np.maximum(result.x, 0)
    return weights / weights.sum()


bank = pl.read_parquet(inp / cfg["prediction_bank"])
groups_source = pl.read_parquet(inp / cfg["group_source"])
metric_rows, fit_rows, prediction_frames, bootstrap_rows = [], [], [], []
rng = np.random.default_rng(cfg["bootstrap_seed"])
for endpoint, settings in cfg["endpoint_settings"].items():
    local = bank.filter(pl.col("endpoint") == endpoint)
    identity = local.select("compound_id", "fold", "y_true").unique().sort("compound_id")
    models = local["model"].unique(maintain_order=True).to_list()
    if models[0] != "current_fixed_ratio":
        raise ValueError("The first prediction-bank model must be the anchor")
    columns = []
    for model in models:
        values = local.filter(pl.col("model") == model).select("compound_id", "prediction")
        joined = identity.select("compound_id").join(
            values, on="compound_id", validate="1:1", maintain_order="left"
        )
        columns.append(joined["prediction"].to_numpy())
    x = np.column_stack(columns)
    folds = identity["fold"].to_numpy().astype(int)
    y = identity["y_true"].to_numpy().astype(int)
    nonanchor_mean = x[:, 1:].mean(axis=1)
    predictions = {"current_fixed_ratio": x[:, 0]}
    for method in ["anchor_mean", "anchor_pairwise_simplex"]:
        prediction = np.empty(len(y), dtype=float)
        for held in np.unique(folds):
            train, valid = folds != held, folds == held
            if method == "anchor_mean":
                candidates = [
                    (1 - alpha) * x[:, 0] + alpha * nonanchor_mean for alpha in cfg["mean_alpha"]
                ]
                values = [score(y[train], p[train], folds[train], settings) for p in candidates]
                choice = int(np.argmax(values))
                prediction[valid] = candidates[choice][valid]
                selected = cfg["mean_alpha"][choice]
                weights = None
            else:
                candidate_weights = [
                    fit_anchor_simplex(x[train], y[train], value) for value in cfg["anchor_lambda"]
                ]
                values = [
                    score(y[train], (x @ weights)[train], folds[train], settings)
                    for weights in candidate_weights
                ]
                choice = int(np.argmax(values))
                weights = candidate_weights[choice]
                prediction[valid] = x[valid] @ weights
                selected = cfg["anchor_lambda"][choice]
            fit_rows.append(
                {
                    "endpoint": endpoint,
                    "method": method,
                    "held_fold": int(held),
                    "selected": selected,
                    "inner_scores": json.dumps(values),
                    "weights": None
                    if weights is None
                    else json.dumps(dict(zip(models, weights.tolist(), strict=True))),
                }
            )
        predictions[method] = prediction
    calls_by_method = {}
    for method, prediction in predictions.items():
        calls = fixed_calls(prediction, folds, settings["call_rate"])
        calls_by_method[method] = calls
        metric_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls, settings["target_prior"]),
            }
        )
        prediction_frames.append(
            identity.with_columns(
                pl.lit(method).alias("method"),
                pl.Series("prediction", prediction),
                pl.Series("call", calls),
            )
        )
    group_table = groups_source.filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["group_method"])
    ).select("compound_id", pl.col("group").cast(pl.String))
    groups = (
        identity.select("compound_id")
        .join(group_table, on="compound_id", validate="1:1", maintain_order="left")["group"]
        .to_numpy()
    )
    group_values = np.unique(groups)
    group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
    for method in ["anchor_mean", "anchor_pairwise_simplex"]:
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows = np.concatenate([group_rows[group] for group in sampled])
            deltas.append(
                expected_mcc(y[rows], calls_by_method[method][rows], settings["target_prior"])
                - expected_mcc(
                    y[rows],
                    calls_by_method["current_fixed_ratio"][rows],
                    settings["target_prior"],
                )
            )
        delta = np.asarray(deltas)
        bootstrap_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "mean_MCC_delta": float(delta.mean()),
                "ci025": float(np.quantile(delta, 0.025)),
                "ci975": float(np.quantile(delta, 0.975)),
                "probability_positive": float((delta > 0).mean()),
            }
        )

metrics = pl.DataFrame(metric_rows).sort(["endpoint", "MCC_target_prior"], descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(fit_rows).write_parquet(out / "fits.parquet")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.concat(prediction_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(json.dumps({"protocol": cfg["protocol"]}, indent=2) + "\n")
print(metrics, flush=True)
