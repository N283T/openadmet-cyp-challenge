"""Compare old and PXR-extended nested Caruana at fixed A4 call rate."""

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


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values))
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask]) / (mask.sum() + 1.0)
    return result


def calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        result[positions[np.argsort(-values[positions], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    tpr = float(selected[positive].mean())
    fpr = float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


oof_tables = {"old": read("old_oof"), "new": read("new_oof")}
test_tables = {"old": read("old_test"), "new": read("new_test")}
rows = []
test_ranks = {}
for version in ["old", "new"]:
    for method in cfg["methods"]:
        oof = (
            oof_tables[version]
            .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == method))
            .sort("compound_id")
        )
        folds = oof["fold"].to_numpy()
        y = oof["y_true"].to_numpy().astype(int)
        prediction = fold_rank(oof["prediction"].to_numpy(), folds)
        selected = calls(prediction, folds)
        test = (
            test_tables[version]
            .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == method))
            .sort("compound_id")
        )
        test_rank = rankdata(test["prediction"].to_numpy()) / (test.height + 1.0)
        test_ranks[(version, method)] = (test["compound_id"].to_list(), test_rank)
        rows.append(
            {
                "version": version,
                "method": method,
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "MCC_observed": float(matthews_corrcoef(y, selected)),
                "MCC_target_prior": expected_mcc(y, selected),
            }
        )

comparisons = []
for method in cfg["methods"]:
    old_ids, old_rank = test_ranks[("old", method)]
    new_ids, new_rank = test_ranks[("new", method)]
    if old_ids != new_ids:
        raise ValueError(f"Test IDs differ for {method}")
    old_calls = set(np.argsort(-old_rank, kind="stable")[: cfg["test_calls"]])
    new_calls = set(np.argsort(-new_rank, kind="stable")[: cfg["test_calls"]])
    comparisons.append(
        {
            "method": method,
            "test_spearman_old_new": float(spearmanr(old_rank, new_rank).statistic),
            "test_swaps_each_way": len(old_calls - new_calls),
        }
    )

metrics = pl.DataFrame(rows).sort("method", "version")
metrics.write_csv(out / "fixed-rate-metrics.csv")
pl.DataFrame(comparisons).write_csv(out / "test-comparison.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "scope": "Matched fixed-rate evaluation of nested recipe outputs",
            "submission": "None",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
print(pl.DataFrame(comparisons), flush=True)
