"""Measure prediction dependence and leave-one-out effects in selected D6 components."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    tpr = float(selected[positive].mean())
    fpr = float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


reference = (
    read("reference_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["reference_method"]))
    .sort("compound_id")
)
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)


def aligned(source: str, method: str, column: str) -> np.ndarray:
    table = (
        read(source)
        .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
        .select("compound_id", pl.col(column).alias("prediction"))
    )
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if joined.height != ids.height or not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete alignment for {source}/{method}")
    return joined["prediction"].to_numpy()


specs = {
    "lgbm::morgan_plain": ("lgbm_oof", "morgan_plain", "prediction"),
    "lgbm::mordred_plain": ("lgbm_oof", "mordred_plain", "prediction"),
    "lgbm::fp_mordred_plain": ("lgbm_oof", "fp_mordred_plain", "prediction"),
    "lgbm::external_class_plain": ("lgbm_oof", "external_class_plain", "prediction"),
    "lgbm::external_mixed_plain": ("lgbm_oof", "external_mixed_plain", "prediction"),
    "pfn::external_mixed_pca256": ("pfn_oof", "external_mixed_pca256", "probability"),
    "graph::tdi_only": ("graph_oof", "tdi_only", "probability"),
}
raw_core = {name: aligned(*specs[name]) for name in cfg["selected_core"]}
components = {name: fold_rank(value, folds) for name, value in raw_core.items()}
fixed_table = read("fixed_oof").filter(pl.col("endpoint") == endpoint)
for name, method in [
    ("morgan::logistic", "logistic_c0.05"),
    ("morgan::rbf", "rbf_svm_c1"),
    ("morgan::complement_nb", "complement_nb_a1"),
]:
    table = fixed_table.filter(pl.col("method") == method).select(
        "compound_id", "prediction"
    )
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    components[name] = fold_rank(joined["prediction"].to_numpy(), folds)


def blend(names: list[str]) -> np.ndarray:
    core_names = [name for name in names if name in raw_core]
    core = fold_rank(np.mean([raw_core[name] for name in core_names], axis=0), folds)
    weights = {
        "core": 4 / 9,
        "morgan::logistic": 2 / 9,
        "morgan::rbf": 1 / 6,
        "morgan::complement_nb": 1 / 6,
    }
    available = {"core": core} | {
        name: components[name] for name in weights if name != "core" and name in names
    }
    denominator = sum(weights[name] for name in available)
    return sum(weights[name] * value for name, value in available.items()) / denominator


all_names = list(components)
full = blend(all_names)
full_calls = calls(full, folds)
metric_rows = []
for name, prediction in components.items():
    selected = calls(prediction, folds)
    metric_rows.append(
        {
            "component": name,
            "contains_morgan": "morgan" in name or name == "lgbm::fp_mordred_plain",
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_target_prior": expected_mcc(y, selected),
            "spearman_with_full": float(spearmanr(prediction, full).statistic),
        }
    )

pair_rows = []
for index, left in enumerate(all_names):
    for right in all_names[index + 1 :]:
        left_calls, right_calls = calls(components[left], folds), calls(components[right], folds)
        pair_rows.append(
            {
                "left": left,
                "right": right,
                "spearman": float(spearmanr(components[left], components[right]).statistic),
                "call_jaccard": float(
                    np.logical_and(left_calls, right_calls).sum()
                    / np.logical_or(left_calls, right_calls).sum()
                ),
            }
        )

ablation_rows = []
for removed in all_names:
    prediction = blend([name for name in all_names if name != removed])
    selected = calls(prediction, folds)
    ablation_rows.append(
        {
            "removed": removed,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_target_prior": expected_mcc(y, selected),
            "delta_AUROC": float(roc_auc_score(y, prediction) - roc_auc_score(y, full)),
            "delta_AP": float(
                average_precision_score(y, prediction) - average_precision_score(y, full)
            ),
            "delta_MCC": expected_mcc(y, selected) - expected_mcc(y, full_calls),
        }
    )

full_row = {
    "AUROC": float(roc_auc_score(y, full)),
    "AP": float(average_precision_score(y, full)),
    "MCC_target_prior": expected_mcc(y, full_calls),
}
pl.DataFrame(metric_rows).write_csv(out / "component-metrics.csv")
pl.DataFrame(pair_rows).write_csv(out / "pairwise-correlations.csv")
pl.DataFrame(ablation_rows).write_csv(out / "leave-one-out.csv")
(out / "summary.json").write_text(
    json.dumps({"fixed_full_oof": full_row, "selection": cfg["selection"]}, indent=2) + "\n"
)
print(json.dumps(full_row), flush=True)
