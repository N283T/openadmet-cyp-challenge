"""Audit BioTransformer product-activity ranks against the current live A4 rank."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
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


def calls(values: np.ndarray, folds: np.ndarray | None = None) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    if folds is None:
        result[np.argsort(-values, kind="stable")[: cfg["positive_calls"]]] = True
        return result
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = int(round(cfg["call_rate"] * len(rows)))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


caruana = read("caruana_oof").sort("compound_id")
folds = caruana["fold"].to_numpy().astype(int)
y = caruana["y_true"].to_numpy().astype(int)
groups = caruana.select("compound_id").join(
    read("activity_oof")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["candidate_methods"][0])
    )
    .select("compound_id", "group"),
    on="compound_id",
    validate="1:1",
)["group"].cast(pl.String).to_numpy()
pxr = (
    read("pxr_oof")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["pxr_method"])
    )
    .select("compound_id", "prediction")
)
pxr = caruana.select("compound_id").join(pxr, on="compound_id", validate="1:1")
current_rank = (1 - cfg["pxr_weight"]) * fold_rank(
    caruana["caruana_auc"].to_numpy(), folds
) + cfg["pxr_weight"] * fold_rank(pxr["prediction"].to_numpy(), folds)

current_test = (
    read("current_test")
    .filter(pl.col("method") == cfg["current_test_method"])
    .select("compound_id", "prediction")
    .sort("compound_id")
)
current_test_rank = unit_rank(current_test["prediction"].to_numpy())
current_calls = calls(current_rank, folds)
current_test_calls = calls(current_test_rank)

metric_rows = [
    {
        "method": "current",
        "weight": 0.0,
        "AUROC": float(roc_auc_score(y, current_rank)),
        "AP": float(average_precision_score(y, current_rank)),
        "MCC_observed": float(matthews_corrcoef(y, current_calls)),
        "MCC_target_prior": expected_mcc(y, current_calls),
        "oof_spearman_to_current": 1.0,
        "test_spearman_to_current": 1.0,
        "test_swaps": 0,
    }
]
bootstrap_rows = []
test_rows = []
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
for method in cfg["candidate_methods"]:
    candidate = (
        read("activity_oof")
        .filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == "umap")
            & (pl.col("method") == method)
        )
        .select("compound_id", "prediction")
    )
    candidate = caruana.select("compound_id").join(candidate, on="compound_id", validate="1:1")
    candidate_rank = fold_rank(candidate["prediction"].to_numpy(), folds)
    candidate_test = (
        read("activity_test")
        .filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == "umap")
            & (pl.col("method") == method)
        )
        .select("compound_id", "prediction")
    )
    candidate_test = current_test.select("compound_id").join(
        candidate_test, on="compound_id", validate="1:1"
    )
    candidate_test_rank = unit_rank(candidate_test["prediction"].to_numpy())
    for weight in cfg["weights"]:
        blend = (1 - weight) * current_rank + weight * candidate_rank
        test_blend = (1 - weight) * current_test_rank + weight * candidate_test_rank
        selected = calls(blend, folds)
        test_selected = calls(test_blend)
        name = f"current+{weight:g}*{method}"
        metric_rows.append(
            {
                "method": method,
                "weight": weight,
                "AUROC": float(roc_auc_score(y, blend)),
                "AP": float(average_precision_score(y, blend)),
                "MCC_observed": float(matthews_corrcoef(y, selected)),
                "MCC_target_prior": expected_mcc(y, selected),
                "oof_spearman_to_current": float(
                    spearmanr(current_rank, candidate_rank).statistic
                ),
                "test_spearman_to_current": float(
                    spearmanr(current_test_rank, candidate_test_rank).statistic
                ),
                "test_swaps": int((current_test_calls & ~test_selected).sum()),
            }
        )
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows = np.concatenate([group_rows[group] for group in sampled])
            deltas.append(
                expected_mcc(y[rows], selected[rows])
                - expected_mcc(y[rows], current_calls[rows])
            )
        bootstrap_rows.append(
            {
                "name": name,
                "delta": expected_mcc(y, selected) - expected_mcc(y, current_calls),
                "p025": float(np.quantile(deltas, 0.025)),
                "p975": float(np.quantile(deltas, 0.975)),
                "improvement_probability": float(np.mean(np.asarray(deltas) > 0)),
            }
        )
        for compound_id, prediction, call in zip(
            current_test["compound_id"], test_blend, test_selected, strict=True
        ):
            test_rows.append(
                {
                    "compound_id": compound_id,
                    "endpoint": endpoint,
                    "method": name,
                    "prediction": float(prediction),
                    "call": bool(call),
                    "fold": -1,
                }
            )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
bootstrap = pl.DataFrame(bootstrap_rows).sort("delta", descending=True)
metrics.write_csv(out / "metrics.csv")
bootstrap.write_csv(out / "conditional-bootstrap.csv")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "best": metrics.row(0, named=True),
            "scope": "Current-live A4 increment audit; no candidate or submission.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
print(bootstrap, flush=True)
