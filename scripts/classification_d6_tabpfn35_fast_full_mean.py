"""Evaluate a matched fast/full TabPFN 3.5 mean for CYP2D6."""

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


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return calls


def top_calls(values: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    calls[np.argsort(-values, kind="stable")[: cfg["test_calls"]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


base = read("base_oof").filter(pl.col("model") == "deployed_mean_block").sort("compound_id")
ids = base.select("compound_id")
folds = base["fold"].to_numpy().astype(int)
groups = base["group"].to_numpy()
y = base["y_true"].to_numpy().astype(int)
base_prediction = fold_rank(base["prediction"].to_numpy(), folds)
base_calls = fixed_calls(base_prediction, folds)


def aligned(source: str, ids_frame: pl.DataFrame, test: bool = False) -> np.ndarray:
    table = read(source).filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == cfg["method"])
    ).select("compound_id", "prediction")
    joined = ids_frame.join(table, on="compound_id", validate="1:1", maintain_order="left")
    values = joined["prediction"].to_numpy()
    return unit_rank(values) if test else fold_rank(values, folds)


fast = aligned("fast_oof", ids)
full = aligned("full_oof", ids)
mean_prediction = (fast + full) / 2

deployed = read("deployed_test").sort("compound_id")
test_ids = deployed.select("compound_id")
base_test = unit_rank(deployed["d6_score"].to_numpy())
base_test_calls = top_calls(base_test)
fast_test = aligned("fast_test", test_ids, test=True)
full_test = aligned("full_test", test_ids, test=True)
mean_test = (fast_test + full_test) / 2

metric_rows = []
test_rows = []
bootstrap_rows = []
models = [("deployed_mean_block", base_prediction, base_test)]
for weight in cfg["weights"]:
    models.append(
        (
            f"tabpfn35_fast_full_mean::w{weight:g}",
            (1 - weight) * base_prediction + weight * mean_prediction,
            (1 - weight) * base_test + weight * mean_test,
        )
    )

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
for name, prediction, test_prediction in models:
    calls = fixed_calls(prediction, folds)
    test_calls = top_calls(test_prediction)
    metric_rows.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "swaps_from_base": int((calls & ~base_calls).sum()),
        }
    )
    test_rows.append(
        {
            "model": name,
            "positive_calls": int(test_calls.sum()),
            "swaps_from_deployed": int((test_calls & ~base_test_calls).sum()),
        }
    )
    if name == "deployed_mean_block":
        continue
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], calls[rows]) - expected_mcc(y[rows], base_calls[rows])
        )
    delta = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "model": name,
            "mean_MCC_delta": float(delta.mean()),
            "ci025": float(np.quantile(delta, 0.025)),
            "ci975": float(np.quantile(delta, 0.975)),
            "probability_positive": float((delta > 0).mean()),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(test_rows).write_csv(out / "test-swaps.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
