"""Audit TabPFN 3.5-fast additions to the scored A4 equal-seven rank."""

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


base = read("base_oof").filter(pl.col("model") == cfg["base_model"]).sort("compound_id")
ids = base.select("compound_id")
folds = base["fold"].to_numpy().astype(int)
y = base["y_true"].to_numpy().astype(int)
groups = (
    ids.join(
        read("assignments")
        .filter(pl.col("strategy") == cfg["split_id"])
        .select("compound_id", "group"),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )["group"]
    .to_numpy()
)
base_prediction = fold_rank(base["prediction"].to_numpy(), folds)
base_calls = fixed_calls(base_prediction, folds)

new = (
    read("new_oof")
    .filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == cfg["new_method"])
    )
    .select("compound_id", "prediction")
)
new_prediction = fold_rank(
    ids.join(new, on="compound_id", validate="1:1", maintain_order="left")[
        "prediction"
    ].to_numpy(),
    folds,
)

base_test = read("base_test").filter(pl.col("model") == cfg["base_model"]).sort("compound_id")
test_ids = base_test.select("compound_id")
base_test_prediction = unit_rank(base_test["prediction"].to_numpy())
base_test_calls = top_calls(base_test_prediction)
new_test = (
    read("new_test")
    .filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == cfg["new_method"])
    )
    .select("compound_id", "prediction")
)
new_test_prediction = unit_rank(
    test_ids.join(new_test, on="compound_id", validate="1:1", maintain_order="left")[
        "prediction"
    ].to_numpy()
)

models = [
    ("scored_equal7", base_prediction, base_test_prediction),
    (cfg["new_label"], new_prediction, new_test_prediction),
]
for weight in cfg["weights"]:
    models.append(
        (
            f"equal7_{cfg['new_label']}::w{weight:g}",
            (1 - weight) * base_prediction + weight * new_prediction,
            (1 - weight) * base_test_prediction + weight * new_test_prediction,
        )
    )

metric_rows = []
test_rows = []
bootstrap_rows = []
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
            "swaps_from_equal7": int((calls & ~base_calls).sum()),
        }
    )
    test_rows.append(
        {
            "model": name,
            "positive_calls": int(test_calls.sum()),
            "swaps_from_scored_equal7": int((test_calls & ~base_test_calls).sum()),
        }
    )
    if name in {"scored_equal7", cfg["new_label"]}:
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
