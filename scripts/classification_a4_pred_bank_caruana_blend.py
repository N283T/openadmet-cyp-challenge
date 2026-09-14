"""Screen conservative rank blends of auxiliary A4 models with live Caruana."""

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


def rank_by_fold(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    ranked = np.empty(len(values))
    for fold in np.unique(folds):
        mask = folds == fold
        ranked[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return ranked


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = int(round(cfg["call_rate"] * len(positions)))
        result[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
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


base_oof = read("caruana_oof").sort("compound_id")
base_test = read("caruana_test").sort("compound_id")
new_oof = read("new_oof").filter(pl.col("strategy") == "umap")
new_test = read("new_test").filter(pl.col("strategy") == "umap")
folds = base_oof["fold"].to_numpy()
y = base_oof["y_true"].to_numpy().astype(int)
base_oof_rank = rank_by_fold(base_oof[cfg["base_method"]].to_numpy(), folds)
base_test_rank = unit_rank(base_test[cfg["base_method"]].to_numpy())
base_test_calls = set(
    base_test.with_columns(pl.Series("score", base_test_rank))
    .sort("score", descending=True)
    .head(round(cfg["call_rate"] * base_test.height))["compound_id"]
)

rows = []
prediction_rows = []
baseline_calls = calls(base_oof_rank, folds)
rows.append(
    {
        "method": cfg["base_method"],
        "weight": 0.0,
        "AUROC": float(roc_auc_score(y, base_oof_rank)),
        "AP": float(average_precision_score(y, base_oof_rank)),
        "MCC_observed": float(matthews_corrcoef(y, baseline_calls)),
        "MCC_target_prior": expected_mcc(y, baseline_calls),
        "oof_spearman_to_base": 1.0,
        "test_spearman_to_base": 1.0,
        "test_call_swaps": 0,
    }
)
for method in cfg["candidate_methods"]:
    candidate_oof = (
        new_oof.filter(pl.col("method") == method)
        .select("compound_id", "prediction")
        .join(base_oof.select("compound_id"), on="compound_id", maintain_order="right")
    )
    candidate_test = (
        new_test.filter(pl.col("method") == method)
        .select("compound_id", "prediction")
        .join(base_test.select("compound_id"), on="compound_id", maintain_order="right")
    )
    candidate_oof_rank = rank_by_fold(candidate_oof["prediction"].to_numpy(), folds)
    candidate_test_rank = unit_rank(candidate_test["prediction"].to_numpy())
    for weight in cfg["weights"]:
        blended_oof = (1 - weight) * base_oof_rank + weight * candidate_oof_rank
        blended_test = (1 - weight) * base_test_rank + weight * candidate_test_rank
        selected = calls(blended_oof, folds)
        test_selected = set(
            base_test.with_columns(pl.Series("score", blended_test))
            .sort("score", descending=True)
            .head(round(cfg["call_rate"] * base_test.height))["compound_id"]
        )
        name = f"{cfg['base_method']}+{weight:g}*{method}"
        rows.append(
            {
                "method": method,
                "weight": weight,
                "AUROC": float(roc_auc_score(y, blended_oof)),
                "AP": float(average_precision_score(y, blended_oof)),
                "MCC_observed": float(matthews_corrcoef(y, selected)),
                "MCC_target_prior": expected_mcc(y, selected),
                "oof_spearman_to_base": float(
                    spearmanr(base_oof_rank, candidate_oof_rank).statistic
                ),
                "test_spearman_to_base": float(
                    spearmanr(base_test_rank, candidate_test_rank).statistic
                ),
                "test_call_swaps": len(base_test_calls - test_selected),
            }
        )
        for index, compound_id in enumerate(base_test["compound_id"]):
            prediction_rows.append(
                {
                    "compound_id": compound_id,
                    "method": name,
                    "prediction": float(blended_test[index]),
                }
            )

screen = pl.DataFrame(rows).sort("MCC_target_prior", descending=True)
screen.write_csv(out / "screen.csv")
pl.DataFrame(prediction_rows).write_parquet(out / "test-blends.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "baseline": cfg["base_method"],
            "selection_warning": (
                "Descriptive same-OOF screen only; gains require cross-fit or bootstrap support "
                "before deployment consideration."
            ),
            "best": screen.row(0, named=True),
        },
        indent=2,
    )
    + "\n"
)
print(screen, flush=True)
