"""Test whether the dual-split A4 block complements the live Caruana rank."""

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
endpoint = "CYP3A4_is_TDI"


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


def top_calls(scores: np.ndarray, folds: np.ndarray | None = None) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    if folds is None:
        calls[np.argsort(-scores, kind="stable")[: cfg["positive_calls"]]] = True
        return calls
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = int(round(cfg["call_rate"] * len(rows)))
        order = np.argsort(-scores[rows], kind="stable")
        calls[rows[order[:count]]] = True
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


caruana = (
    read("caruana_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["caruana_method"]))
    .sort("compound_id")
)
block = (
    read("block_oof")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["block_method"])
    )
    .sort("compound_id")
)
keys = ["compound_id", "fold", "y_true"]
if not caruana.select(keys).equals(block.select(keys)):
    raise ValueError("Caruana and block OOF keys differ")
y = block["y_true"].to_numpy().astype(int)
folds = block["fold"].to_numpy().astype(int)
groups = block["group"].cast(pl.String).to_numpy()
caruana_rank = fold_rank(caruana["prediction"].to_numpy(), folds)
block_rank = fold_rank(block["prediction"].to_numpy(), folds)

caruana_test = read("caruana_test").filter(
    (pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["caruana_method"])
)
if caruana_test["compound_id"].n_unique() != len(caruana_test):
    caruana_test = caruana_test.group_by("compound_id").agg(pl.col("prediction").mean())
caruana_test = caruana_test.sort("compound_id")
block_test = (
    read("block_test")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["block_method"])
    )
    .sort("compound_id")
)
if caruana_test["compound_id"].to_list() != block_test["compound_id"].to_list():
    raise ValueError("Caruana and block test keys differ")
caruana_test_rank = unit_rank(caruana_test["prediction"].to_numpy())
block_test_rank = unit_rank(block_test["prediction"].to_numpy())

predictions = {}
test_predictions = {}
metric_rows = []
for method, weights in cfg["weights"].items():
    predictions[method] = (weights[0] * caruana_rank + weights[1] * block_rank) / sum(weights)
    test_predictions[method] = (
        weights[0] * caruana_test_rank + weights[1] * block_test_rank
    ) / sum(weights)
    calls = top_calls(predictions[method], folds)
    metric_rows.append(
        {
            "method": method,
            "caruana_weight": weights[0],
            "block_weight": weights[1],
            "AUROC": float(roc_auc_score(y, predictions[method])),
            "AP": float(average_precision_score(y, predictions[method])),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
        }
    )

baseline_calls = top_calls(predictions["caruana"], folds)
baseline_test_calls = top_calls(test_predictions["caruana"])
rng = np.random.default_rng(cfg["bootstrap_seed"])
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
bootstrap_rows = []
test_rows = []
for method, prediction in predictions.items():
    calls = top_calls(prediction, folds)
    test_calls = top_calls(test_predictions[method])
    if method != "caruana":
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
                "delta_vs_caruana": expected_mcc(y, calls) - expected_mcc(y, baseline_calls),
                "bootstrap_p025": float(np.quantile(deltas, 0.025)),
                "bootstrap_p975": float(np.quantile(deltas, 0.975)),
                "test_removed": int((baseline_test_calls & ~test_calls).sum()),
                "test_added": int((~baseline_test_calls & test_calls).sum()),
            }
        )
    for row, value in enumerate(test_predictions[method]):
        test_rows.append(
            {
                "compound_id": block_test["compound_id"][row],
                "endpoint": endpoint,
                "method": method,
                "prediction": float(value),
                "call": bool(test_calls[row]),
                "fold": -1,
            }
        )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).sort("delta_vs_caruana", descending=True).write_csv(
    out / "conditional-bootstrap.csv"
)
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "scope": "Post-screen A4 complement analysis; no candidate or submission.",
            "caruana_method": cfg["caruana_method"],
            "block_method": cfg["block_method"],
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
print(pl.DataFrame(bootstrap_rows).write_csv(), flush=True)
