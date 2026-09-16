"""Distill the currently scored D6 and A4 rank ensembles."""

import itertools
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.optimize import minimize
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


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


def calls_by_fold(values: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    result[np.argsort(-values, kind="stable")[:count]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def simplex(features: np.ndarray, target: np.ndarray) -> np.ndarray:
    if features.shape[1] == 1:
        return np.ones(1)
    initial = np.full(features.shape[1], 1 / features.shape[1])
    result = minimize(
        lambda weight: float(np.mean((features @ weight - target) ** 2)),
        initial,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * features.shape[1],
        constraints={"type": "eq", "fun": lambda weight: weight.sum() - 1.0},
        options={"ftol": 1e-14, "maxiter": 2000},
    )
    if not result.success:
        raise RuntimeError(result.message)
    weight = np.clip(result.x, 0, 1)
    return weight / weight.sum()


def evaluate(
    endpoint: str,
    names: list[str],
    oof_matrix: np.ndarray,
    test_matrix: np.ndarray,
    target_oof: np.ndarray,
    target_test: np.ndarray,
    y: np.ndarray,
    folds: np.ndarray,
    rate: float,
    test_count: int,
    subsets: list[tuple[int, ...]],
    target_prior: float | None = None,
) -> tuple[pl.DataFrame, dict[int, np.ndarray]]:
    target_oof_calls = calls_by_fold(target_oof, folds, rate)
    target_test_calls = top_calls(target_test, test_count)
    rows = []
    test_predictions = {}
    for subset in subsets:
        columns = list(subset)
        weight = simplex(oof_matrix[:, columns], target_oof)
        prediction = oof_matrix[:, columns] @ weight
        deploy = test_matrix[:, columns] @ weight
        selected = calls_by_fold(prediction, folds, rate)
        deploy_selected = top_calls(deploy, test_count)
        selected_names = [names[column] for column in columns]
        size = len(columns)
        row = {
            "endpoint": endpoint,
            "requested_size": size,
            "effective_size": int((weight > cfg["effective_weight_min"]).sum()),
            "members": " + ".join(selected_names),
            "weights": json.dumps(dict(zip(selected_names, weight.tolist(), strict=True))),
            "oof_rmse_to_current": float(np.sqrt(np.mean((prediction - target_oof) ** 2))),
            "oof_spearman_to_current": float(spearmanr(prediction, target_oof).statistic),
            "oof_AUROC": float(roc_auc_score(y, prediction)),
            "oof_AP": float(average_precision_score(y, prediction)),
            "oof_MCC": float(matthews_corrcoef(y, selected)),
            "oof_call_swaps_each_way": int((selected & ~target_oof_calls).sum()),
            "test_spearman_to_current": float(spearmanr(deploy, target_test).statistic),
            "test_call_swaps_each_way": int((deploy_selected & ~target_test_calls).sum()),
        }
        if target_prior is not None:
            row["oof_target_prior_MCC"] = expected_mcc(y, selected, target_prior)
        rows.append(row)
        test_predictions[size] = deploy
    return pl.DataFrame(rows), test_predictions


# A4: six effective Caruana members plus the PXR/Admetica auxiliary member.
a4_oof = read("a4_member_oof").sort("compound_id")
a4_test = read("a4_member_test").sort("compound_id")
a4_aux_oof = (
    read("a4_aux_oof")
    .filter(
        (pl.col("endpoint") == "CYP3A4_is_TDI")
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["a4_aux_method"])
    )
    .select("compound_id", "prediction")
)
a4_aux_oof = a4_oof.select("compound_id").join(
    a4_aux_oof, on="compound_id", validate="1:1"
)
a4_aux_test = (
    read("a4_aux_test")
    .filter(
        (pl.col("endpoint") == "CYP3A4_is_TDI")
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["a4_aux_method"])
    )
    .select("compound_id", "prediction")
)
a4_aux_test = a4_test.select("compound_id").join(
    a4_aux_test, on="compound_id", validate="1:1"
)
a4_folds = a4_oof["fold"].to_numpy().astype(int)
a4_y = a4_oof["y_true"].to_numpy().astype(int)
a4_members = [
    column
    for column in a4_oof.columns
    if column not in {"compound_id", "fold", "y_true", "caruana_auc"}
]
a4_names = [*a4_members, "pxr_admetica"]
a4_oof_matrix = np.column_stack(
    [fold_rank(a4_oof[name].to_numpy(), a4_folds) for name in a4_members]
    + [fold_rank(a4_aux_oof["prediction"].to_numpy(), a4_folds)]
)
a4_test_matrix = np.column_stack(
    [unit_rank(a4_test[name].to_numpy()) for name in a4_members]
    + [unit_rank(a4_aux_test["prediction"].to_numpy())]
)
a4_target_oof = cfg["a4_caruana_weight"] * fold_rank(
    a4_oof["caruana_auc"].to_numpy(), a4_folds
) + (1 - cfg["a4_caruana_weight"]) * fold_rank(
    a4_aux_oof["prediction"].to_numpy(), a4_folds
)
a4_target_test = cfg["a4_caruana_weight"] * unit_rank(
    a4_test["caruana_auc"].to_numpy()
) + (1 - cfg["a4_caruana_weight"]) * unit_rank(a4_aux_test["prediction"].to_numpy())
a4_subsets = [
    subset
    for size in range(1, len(a4_names) + 1)
    for subset in itertools.combinations(range(len(a4_names)), size)
]
a4_all, _ = evaluate(
    "CYP3A4_is_TDI",
    a4_names,
    a4_oof_matrix,
    a4_test_matrix,
    a4_target_oof,
    a4_target_test,
    a4_y,
    a4_folds,
    cfg["a4_rate"],
    cfg["a4_calls"],
    a4_subsets,
)

# D6: expand the old eleven-member mean and append Morgan logistic and ComplementNB.
d6_reference_oof = (
    read("d6_reference_oof")
    .filter(
        (pl.col("endpoint") == "CYP2D6_is_TDI")
        & (pl.col("method") == cfg["d6_reference_method"])
    )
    .sort("compound_id")
)
d6_reference_test = (
    read("d6_reference_test")
    .filter(
        (pl.col("endpoint") == "CYP2D6_is_TDI")
        & (pl.col("method") == cfg["d6_reference_method"])
    )
    .sort("compound_id")
)
d6_folds = d6_reference_oof["fold"].to_numpy().astype(int)
d6_y = d6_reference_oof["y_true"].to_numpy().astype(int)


def aligned_oof(source: str, method: str, prediction: str) -> np.ndarray:
    table = (
        read(source)
        .filter(
            (pl.col("endpoint") == "CYP2D6_is_TDI")
            & (pl.col("method") == method)
        )
        .select("compound_id", pl.col(prediction).alias("prediction"))
    )
    joined = d6_reference_oof.select("compound_id").join(
        table, on="compound_id", validate="1:1"
    )
    return fold_rank(joined["prediction"].to_numpy(), d6_folds)


def aligned_test(source: str, method: str, prediction: str) -> np.ndarray:
    table = (
        read(source)
        .filter(
            (pl.col("endpoint") == "CYP2D6_is_TDI")
            & (pl.col("method") == method)
        )
        .group_by("compound_id")
        .agg(pl.col(prediction).mean().alias("prediction"))
    )
    joined = d6_reference_test.select("compound_id").join(
        table, on="compound_id", validate="1:1"
    )
    return unit_rank(joined["prediction"].to_numpy())


d6_specs = [
    ("lgbm::morgan_plain", "lgbm_oof", "lgbm_test", "morgan_plain", "prediction"),
    ("lgbm::mordred_plain", "lgbm_oof", "lgbm_test", "mordred_plain", "prediction"),
    ("lgbm::fp_mordred_plain", "lgbm_oof", "lgbm_test", "fp_mordred_plain", "prediction"),
    (
        "lgbm::external_class_plain",
        "lgbm_oof",
        "lgbm_test",
        "external_class_plain",
        "prediction",
    ),
    (
        "lgbm::external_mixed_plain",
        "lgbm_oof",
        "lgbm_test",
        "external_mixed_plain",
        "prediction",
    ),
    (
        "lgbm::external_class_mordred_plain",
        "lgbm_oof",
        "lgbm_test",
        "external_class_mordred_plain",
        "prediction",
    ),
    (
        "pfn::external_class_pca256",
        "pfn_oof",
        "pfn_test",
        "external_class_pca256",
        "probability",
    ),
    (
        "pfn::external_mixed_pca256",
        "pfn_oof",
        "pfn_test",
        "external_mixed_pca256",
        "probability",
    ),
    (
        "pfn::chemeleon_pca256",
        "pfn_oof",
        "pfn_test",
        "chemeleon_pca256",
        "probability",
    ),
    ("graph::tdi_only", "graph_oof", "graph_test", "tdi_only", "probability"),
    ("graph::tdi_aux", "graph_oof", "graph_test", "tdi_aux", "probability"),
]
d6_names = [spec[0] for spec in d6_specs] + ["morgan_logistic", "complement_nb"]
d6_oof_columns = [
    aligned_oof(oof_source, method, prediction)
    for _, oof_source, _, method, prediction in d6_specs
]
d6_test_columns = [
    aligned_test(test_source, method, prediction)
    for _, _, test_source, method, prediction in d6_specs
]
fixed_oof = read("fixed_oof").sort("compound_id")
for method in ["logistic_c0.05", "complement_nb_a1"]:
    table = fixed_oof.filter(pl.col("method") == method)
    if not table["compound_id"].equals(d6_reference_oof["compound_id"]):
        raise ValueError(f"Misaligned fixed-model OOF: {method}")
    d6_oof_columns.append(fold_rank(table["prediction"].to_numpy(), d6_folds))

compounds = read("compounds").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
if not test["compound_id"].equals(d6_reference_test["compound_id"]):
    raise ValueError("D6 test identities differ")
id_to_index = {compound_id: index for index, compound_id in enumerate(compounds["compound_id"])}
test_indices = np.asarray([id_to_index[compound_id] for compound_id in test["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
for key in ["logistic_models", "complement_models"]:
    fold_predictions = []
    for model_name in cfg[key]:
        model = joblib.load(inp / model_name)
        if hasattr(model, "predict_proba"):
            score = model.predict_proba(morgan[test_indices])[:, 1]
        else:
            score = model.decision_function(morgan[test_indices])
        fold_predictions.append(unit_rank(score))
    d6_test_columns.append(unit_rank(np.mean(fold_predictions, axis=0)))

d6_oof_matrix = np.column_stack(d6_oof_columns)
d6_test_matrix = np.column_stack(d6_test_columns)
d6_old_oof = fold_rank(d6_reference_oof["prediction"].to_numpy(), d6_folds)
d6_old_test = unit_rank(d6_reference_test["prediction"].to_numpy())
d6_live_oof = (2 * d6_old_oof + d6_oof_matrix[:, -2]) / 3
d6_live_test = (2 * d6_old_test + d6_test_matrix[:, -2]) / 3
d6_target_oof = (2 * d6_live_oof + d6_oof_matrix[:, -1]) / 3
d6_target_test = (2 * d6_live_test + d6_test_matrix[:, -1]) / 3

# Exhaust sizes one to three, then extend the best triple greedily.
d6_subsets = [
    subset
    for size in range(1, 4)
    for subset in itertools.combinations(range(len(d6_names)), size)
]
d6_small, _ = evaluate(
    "CYP2D6_is_TDI",
    d6_names,
    d6_oof_matrix,
    d6_test_matrix,
    d6_target_oof,
    d6_target_test,
    d6_y,
    d6_folds,
    cfg["d6_rate"],
    cfg["d6_calls"],
    d6_subsets,
    cfg["d6_target_prior"],
)
best_three = d6_small.filter(pl.col("requested_size") == 3).sort(
    "oof_rmse_to_current"
).row(0, named=True)
selected = tuple(d6_names.index(name) for name in best_three["members"].split(" + "))
d6_greedy_tables = [d6_small]
for _size in range(4, len(d6_names) + 1):
    candidates = [
        tuple(sorted((*selected, candidate)))
        for candidate in range(len(d6_names))
        if candidate not in selected
    ]
    table, _ = evaluate(
        "CYP2D6_is_TDI",
        d6_names,
        d6_oof_matrix,
        d6_test_matrix,
        d6_target_oof,
        d6_target_test,
        d6_y,
        d6_folds,
        cfg["d6_rate"],
        cfg["d6_calls"],
        candidates,
        cfg["d6_target_prior"],
    )
    selected_row = table.sort("oof_rmse_to_current").row(0, named=True)
    selected = tuple(d6_names.index(name) for name in selected_row["members"].split(" + "))
    d6_greedy_tables.append(table.filter(pl.col("members") == selected_row["members"]))
d6_all = pl.concat(d6_greedy_tables, how="diagonal_relaxed")

all_results = pl.concat([a4_all, d6_all], how="diagonal_relaxed")
all_results.write_csv(out / "all-subsets.csv")
best = (
    all_results.sort("endpoint", "requested_size", "oof_rmse_to_current")
    .group_by("endpoint", "requested_size", maintain_order=True)
    .first()
    .sort("endpoint", "requested_size")
)
best.write_csv(out / "best-by-size.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "a4_members": a4_names,
            "d6_members": d6_names,
            "a4_search": "exhaustive all subsets",
            "d6_search": "exhaustive through size three, then greedy forward",
            "caveat": "Distillation fidelity is not independent performance validation.",
        },
        indent=2,
    )
    + "\n"
)
print(best, flush=True)
