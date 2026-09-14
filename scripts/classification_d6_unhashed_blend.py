"""Blend the live D6 rank with the dual-split-supported unhashed-radius rank."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = int(round(cfg["call_rate"] * len(positions)))
        order = np.argsort(-scores[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


radius = (
    read("radius_oof")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["radius_method"])
    )
    .sort("compound_id")
)
live = (
    read("live_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["live_method"]))
    .sort("compound_id")
)
if live["compound_id"].to_list() != radius["compound_id"].to_list():
    raise ValueError("Live and radius OOF keys differ")

y = radius["y_true"].to_numpy().astype(int)
folds = radius["fold"].to_numpy().astype(int)
groups = radius["group"].cast(pl.String).to_numpy()
live_rank = fold_rank(live["prediction"].to_numpy(), folds)
radius_rank = fold_rank(radius["prediction"].to_numpy(), folds)

radius_test = (
    read("radius_test")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["radius_method"])
    )
    .sort("compound_id")
)
live_test = (
    read("live_test")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["live_method"]))
    .sort("compound_id")
)
if live_test["compound_id"].to_list() != radius_test["compound_id"].to_list():
    raise ValueError("Live and radius test keys differ")
live_test_rank = unit_rank(live_test["prediction"].to_numpy())
radius_test_rank = unit_rank(radius_test["prediction"].to_numpy())

weights = cfg["weights"]
oof_predictions = {
    name: (values[0] * live_rank + values[1] * radius_rank) / sum(values)
    for name, values in weights.items()
}
test_predictions = {
    name: (values[0] * live_test_rank + values[1] * radius_test_rank) / sum(values)
    for name, values in weights.items()
}

metric_rows = []
calls_by_method = {}
for method, prediction in oof_predictions.items():
    calls = fixed_rate_calls(prediction, folds)
    calls_by_method[method] = calls
    metric_rows.append(
        {
            "method": method,
            "live_weight": weights[method][0],
            "radius_weight": weights[method][1],
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "calls": int(calls.sum()),
        }
    )

baseline_calls = calls_by_method["live"]
bootstrap_rows = []
rng = np.random.default_rng(cfg["bootstrap_seed"])
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
for method, calls in calls_by_method.items():
    if method == "live":
        continue
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], calls[rows]) - expected_mcc(y[rows], baseline_calls[rows])
        )
    bootstrap_rows.append(
        {
            "method": method,
            "delta_vs_live": expected_mcc(y, calls) - expected_mcc(y, baseline_calls),
            "bootstrap_mean": float(np.mean(deltas)),
            "bootstrap_p025": float(np.quantile(deltas, 0.025)),
            "bootstrap_p975": float(np.quantile(deltas, 0.975)),
        }
    )

deployment_rows = []
for method, prediction in test_predictions.items():
    count = cfg["positive_calls"]
    order = np.argsort(-prediction, kind="stable")
    calls = np.zeros(len(prediction), dtype=bool)
    calls[order[:count]] = True
    for row, value in enumerate(prediction):
        deployment_rows.append(
            {
                "compound_id": radius_test["compound_id"][row],
                "endpoint": endpoint,
                "method": method,
                "prediction": float(value),
                "call": bool(calls[row]),
                "fold": -1,
                "split_id": "full_training",
            }
        )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).sort("delta_vs_live", descending=True).write_csv(
    out / "conditional-bootstrap.csv"
)
pl.DataFrame(deployment_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "scope": "Post-radius-screen blend; exploratory and not a submission candidate.",
            "live_method": cfg["live_method"],
            "radius_method": cfg["radius_method"],
            "positive_calls": cfg["positive_calls"],
            "bootstrap": "Conditional UMAP-group resampling of fixed OOF calls.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
print(pl.DataFrame(bootstrap_rows).write_csv(), flush=True)
