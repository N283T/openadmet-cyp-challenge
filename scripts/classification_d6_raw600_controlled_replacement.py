"""Compare PCA, raw, and averaged mixed-embedding blocks without reselection."""

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


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
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


def aligned(source: str, method: str, column: str = "prediction") -> np.ndarray:
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
    (
        "lgbm::external_class_mordred_plain"
    ): ("lgbm_oof", "external_class_mordred_plain", "prediction"),
    "pfn::external_class_pca256": ("pfn_oof", "external_class_pca256", "probability"),
    "pfn::external_mixed_pca256": ("pfn_oof", "external_mixed_pca256", "probability"),
    "pfn::chemeleon_pca256": ("pfn_oof", "chemeleon_pca256", "probability"),
    "graph::tdi_only": ("graph_oof", "tdi_only", "probability"),
    "graph::tdi_aux": ("graph_oof", "tdi_aux", "probability"),
    "pfn::external_mixed_raw600": ("raw600_oof", "external_mixed_raw600", "probability"),
}
atoms = {name: aligned(*source) for name, source in specs.items()}

fixed_table = read("fixed_oof").filter(pl.col("endpoint") == endpoint)
fixed, groups = {}, None
for name, method in [
    ("logistic", "logistic_c0.05"),
    ("rbf", "rbf_svm_c1"),
    ("complement", "complement_nb_a1"),
]:
    table = fixed_table.filter(pl.col("method") == method).select(
        "compound_id", "prediction", "group"
    )
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    fixed[name] = fold_rank(joined["prediction"].to_numpy(), folds)
    if groups is None:
        groups = joined["group"].to_numpy()


def recipe_prediction(members: list[str], recipe: str, mode: str) -> np.ndarray:
    values = []
    for member in members:
        if member == cfg["pca_member"] and mode == "raw600":
            values.append(atoms[cfg["raw_member"]])
        elif member == cfg["pca_member"] and mode == "mean_block":
            values.append((atoms[cfg["pca_member"]] + atoms[cfg["raw_member"]]) / 2)
        else:
            values.append(atoms[member])
    core = fold_rank(np.mean(values, axis=0), folds)
    live = (2 * core + fixed["logistic"]) / 3
    recipes = {
        "core": core,
        "core_logistic": live,
        "tail_rbf": (2 * live + fixed["rbf"]) / 3,
        "tail_complement": (2 * live + fixed["complement"]) / 3,
    }
    recipes["tail_average"] = (recipes["tail_rbf"] + recipes["tail_complement"]) / 2
    return recipes[recipe]


selection = read("selections").filter(pl.col("model") == cfg["selection_model"]).row(
    0, named=True
)
choices = json.loads(selection["selection"])
predictions = {}
for scope in ("crossfit", "full"):
    selected = [
        choice
        for choice in choices
        if (choice["held_fold"] == "full") == (scope == "full")
    ]
    for mode in ("pca256", "raw600", "mean_block"):
        prediction = np.empty(len(y), dtype=float)
        for choice in selected:
            candidate = recipe_prediction(choice["members"], choice["recipe"], mode)
            if scope == "full":
                prediction[:] = candidate
            else:
                mask = folds == int(choice["held_fold"])
                prediction[mask] = candidate[mask]
        predictions[f"{scope}::{mode}"] = prediction

rows, call_map = [], {}
for name, prediction in predictions.items():
    selected = fixed_calls(prediction, folds)
    call_map[name] = selected
    rows.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, selected)),
            "MCC_target_prior": expected_mcc(y, selected),
        }
    )

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
bootstrap_rows = []
for scope in ("crossfit", "full"):
    base = f"{scope}::pca256"
    for mode in ("raw600", "mean_block"):
        name = f"{scope}::{mode}"
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows_sampled = np.concatenate([group_rows[group] for group in sampled])
            deltas.append(
                expected_mcc(y[rows_sampled], call_map[name][rows_sampled])
                - expected_mcc(y[rows_sampled], call_map[base][rows_sampled])
            )
        delta = np.asarray(deltas)
        bootstrap_rows.append(
            {
                "model": name,
                "reference": base,
                "mean_MCC_delta": float(delta.mean()),
                "ci025": float(np.quantile(delta, 0.025)),
                "ci975": float(np.quantile(delta, 0.975)),
                "probability_positive": float((delta > 0).mean()),
            }
        )

pl.DataFrame(rows).write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
(out / "summary.json").write_text(json.dumps({"selection": cfg["selection"]}, indent=2) + "\n")
print(pl.DataFrame(rows), flush=True)
