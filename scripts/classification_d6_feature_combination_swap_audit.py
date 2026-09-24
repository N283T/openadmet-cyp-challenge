"""Audit equal-weight addition and one-block swaps in the deployed D6 recipe."""

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
    return pl.read_parquet(inp / name)


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
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


reference = (
    read("reference-oof.parquet")
    .filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("method") == cfg["reference_method"])
    )
    .sort("compound_id")
)
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)


def aligned(source: str, method: str, column: str = "prediction") -> np.ndarray:
    table = read(source).filter(
        (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == method)
    )
    if source == "candidate-oof.parquet":
        table = table.filter(pl.col("split_id") == cfg["candidate_strategy"])
    table = table.select("compound_id", pl.col(column).alias("prediction"))
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete alignment for {source}/{method}")
    return joined["prediction"].to_numpy()


specs = {
    "lgbm::morgan_plain": ("lgbm-oof.parquet", "morgan_plain", "prediction"),
    "lgbm::mordred_plain": ("lgbm-oof.parquet", "mordred_plain", "prediction"),
    "lgbm::fp_mordred_plain": ("lgbm-oof.parquet", "fp_mordred_plain", "prediction"),
    "lgbm::external_class_plain": (
        "lgbm-oof.parquet",
        "external_class_plain",
        "prediction",
    ),
    "lgbm::external_mixed_plain": (
        "lgbm-oof.parquet",
        "external_mixed_plain",
        "prediction",
    ),
    "lgbm::external_class_mordred_plain": (
        "lgbm-oof.parquet",
        "external_class_mordred_plain",
        "prediction",
    ),
    "pfn::external_class_pca256": (
        "pfn-oof.parquet",
        "external_class_pca256",
        "probability",
    ),
    "pfn::external_mixed_pca256": (
        "pfn-oof.parquet",
        "external_mixed_pca256",
        "probability",
    ),
    "graph::tdi_only": ("graph-oof.parquet", "tdi_only", "probability"),
    "graph::tdi_aux": ("graph-oof.parquet", "tdi_aux", "probability"),
    "pfn::external_mixed_raw600": (
        "raw600-oof.parquet",
        "external_mixed_raw600",
        "probability",
    ),
}
atoms = {name: aligned(*source) for name, source in specs.items()}

fixed_table = read("fixed-oof.parquet").filter(pl.col("endpoint") == cfg["endpoint"])
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


def recipe_prediction(members: list[str], recipe: str, candidate=None, replace=None):
    values = []
    for member in members:
        if member == replace:
            values.append(candidate)
        elif member == cfg["pca_member"]:
            values.append((atoms[member] + atoms[cfg["raw_member"]]) / 2)
        else:
            values.append(atoms[member])
    if candidate is not None and replace is None:
        values.append(candidate)
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


choices = json.loads(
    read("selections.parquet")
    .filter(pl.col("model") == cfg["selection_model"])["selection"][0]
)
full_choice = next(choice for choice in choices if choice["held_fold"] == "full")
crossfit_choices = [choice for choice in choices if choice["held_fold"] != "full"]
control = recipe_prediction(full_choice["members"], full_choice["recipe"])
crossfit_control = np.empty(len(y), dtype=float)
for choice in crossfit_choices:
    mask = folds == int(choice["held_fold"])
    crossfit_control[mask] = recipe_prediction(choice["members"], choice["recipe"])[mask]

candidates_table = read("candidate-oof.parquet").filter(
    (pl.col("endpoint") == cfg["endpoint"])
    & (pl.col("split_id") == cfg["candidate_strategy"])
)
candidate_methods = sorted(candidates_table["method"].unique())
candidate_values = {
    method: aligned("candidate-oof.parquet", method) for method in candidate_methods
}

predictions = {"fixed::control": control, "crossfit::control": crossfit_control}
candidate_corr_rows = []
for method, candidate in candidate_values.items():
    candidate_corr_rows.append(
        {
            "candidate": method,
            "control_spearman": float(spearmanr(candidate, control).statistic),
            **{
                f"corr::{member}": float(spearmanr(candidate, values).statistic)
                for member, values in atoms.items()
                if member in full_choice["members"]
            },
        }
    )
    predictions[f"fixed::add::{method}"] = recipe_prediction(
        full_choice["members"], full_choice["recipe"], candidate=candidate
    )
    for member in full_choice["members"]:
        predictions[f"fixed::replace::{member}::{method}"] = recipe_prediction(
            full_choice["members"],
            full_choice["recipe"],
            candidate=candidate,
            replace=member,
        )
    crossfit = np.empty(len(y), dtype=float)
    for choice in crossfit_choices:
        mask = folds == int(choice["held_fold"])
        crossfit[mask] = recipe_prediction(
            choice["members"], choice["recipe"], candidate=candidate
        )[mask]
    predictions[f"crossfit::add::{method}"] = crossfit

metric_rows, call_map = [], {}
for model, prediction in predictions.items():
    calls = fixed_calls(prediction, folds)
    call_map[model] = calls
    reference_name = "crossfit::control" if model.startswith("crossfit") else "fixed::control"
    reference_calls = fixed_calls(predictions[reference_name], folds)
    metric_rows.append(
        {
            "model": model,
            "reference": reference_name,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "delta_MCC": expected_mcc(y, calls)
            - expected_mcc(y, reference_calls),
            "swapped_calls": int(np.count_nonzero(calls != reference_calls)),
        }
    )

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
bootstrap_indices = [
    np.concatenate([group_rows[group] for group in rng.choice(group_values, len(group_values))])
    for _ in range(cfg["bootstrap_replicates"])
]
bootstrap_rows = []
for row in metric_rows:
    model, reference_name = row["model"], row["reference"]
    if model.endswith("::control"):
        continue
    deltas = np.asarray(
        [
            expected_mcc(y[index], call_map[model][index])
            - expected_mcc(y[index], call_map[reference_name][index])
            for index in bootstrap_indices
        ]
    )
    bootstrap_rows.append(
        {
            "model": model,
            "mean_delta": float(deltas.mean()),
            "ci025": float(np.quantile(deltas, 0.025)),
            "ci975": float(np.quantile(deltas, 0.975)),
            "probability_positive": float((deltas > 0).mean()),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("delta_MCC", descending=True)
bootstrap = pl.DataFrame(bootstrap_rows)
metrics.write_csv(out / "metrics.csv")
bootstrap.write_csv(out / "bootstrap.csv")
pl.DataFrame(candidate_corr_rows).write_csv(out / "candidate-correlations.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "fixed_members": full_choice["members"],
            "best_point_deltas": metrics.head(20).to_dicts(),
        },
        indent=2,
    )
    + "\n"
)
print(metrics.head(30), flush=True)
