"""Blend strict PXR Chemprop readouts on fixed validation folds."""

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


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values) / (len(values) + 1.0)


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


rows = []
test_rows = []
for strategy in cfg["strategies"]:
    base_oof = (
        read("base_oof")
        .filter((pl.col("strategy") == strategy) & (pl.col("method") == cfg["base_method"]))
        .sort("compound_id")
    )
    base_test = (
        read("base_test")
        .filter((pl.col("strategy") == strategy) & (pl.col("method") == cfg["base_method"]))
        .sort("compound_id")
    )
    folds = base_oof["fold"].to_numpy()
    y = base_oof["y_true"].to_numpy().astype(int)
    base_rank = fold_rank(base_oof["prediction"].to_numpy(), folds)
    base_test_rank = unit_rank(base_test["prediction"].to_numpy())
    predictions = {cfg["base_method"]: (base_rank, base_test_rank, 1.0)}
    candidate_ranks = {}
    for method in cfg["candidate_methods"]:
        candidate_oof = (
            read("model_oof")
            .filter((pl.col("strategy") == strategy) & (pl.col("method") == method))
            .sort("compound_id")
        )
        candidate_test = (
            read("model_test")
            .filter((pl.col("strategy") == strategy) & (pl.col("method") == method))
            .sort("compound_id")
        )
        if not base_oof["compound_id"].equals(candidate_oof["compound_id"]):
            raise ValueError(f"OOF IDs differ for {method}")
        candidate_rank = fold_rank(candidate_oof["prediction"].to_numpy(), folds)
        candidate_test_rank = unit_rank(candidate_test["prediction"].to_numpy())
        candidate_ranks[method] = (candidate_rank, candidate_test_rank)
        correlation = float(spearmanr(base_rank, candidate_rank).statistic)
        for weight in cfg["candidate_weights"]:
            name = f"{cfg['base_method']}+{weight:g}*{method}"
            predictions[name] = (
                (1 - weight) * base_rank + weight * candidate_rank,
                (1 - weight) * base_test_rank + weight * candidate_test_rank,
                correlation,
            )
    extra_oof, extra_test = candidate_ranks["extra_trees_leaf2"]
    rbf_oof, rbf_test = candidate_ranks["rbf_svm_c1"]
    for base_weight, extra_weight, rbf_weight in [
        (0.5, 0.25, 0.25),
        (0.6, 0.2, 0.2),
        (0.5, 0.3, 0.2),
    ]:
        name = f"threeway_{base_weight:g}_{extra_weight:g}_{rbf_weight:g}"
        blended_oof = base_weight * base_rank + extra_weight * extra_oof + rbf_weight * rbf_oof
        blended_test = (
            base_weight * base_test_rank + extra_weight * extra_test + rbf_weight * rbf_test
        )
        predictions[name] = (
            blended_oof,
            blended_test,
            float(spearmanr(base_rank, blended_oof).statistic),
        )
    for name, (prediction, test_prediction, correlation) in predictions.items():
        selected = calls(prediction, folds)
        rows.append(
            {
                "strategy": strategy,
                "method": name,
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "MCC_observed": float(matthews_corrcoef(y, selected)),
                "MCC_target_prior": expected_mcc(y, selected),
                "oof_spearman_to_base": correlation,
            }
        )
        for index, compound_id in enumerate(base_test["compound_id"]):
            test_rows.append(
                {
                    "compound_id": compound_id,
                    "strategy": strategy,
                    "method": name,
                    "prediction": float(test_prediction[index]),
                    "fold": -1,
                }
            )

metrics = pl.DataFrame(rows).sort("strategy", "MCC_target_prior", descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "scope": "Small same-OOF model-family blend screen",
            "selection_warning": "Descriptive only; no deployment selection.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
