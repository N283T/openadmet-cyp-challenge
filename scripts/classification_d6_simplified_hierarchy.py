"""Evaluate pruned D6 cores after re-optimizing a small downstream recipe."""

import json
import os
from itertools import combinations
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


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
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
    tpr = float(selected[positive].mean())
    fpr = float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def score_objective(
    y: np.ndarray, prediction: np.ndarray, folds: np.ndarray, objective: str
) -> float:
    if objective == "auc":
        return float(roc_auc_score(y, prediction))
    return expected_mcc(y, fixed_calls(prediction, folds))


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


specs = [
    ("lgbm::morgan_plain", "lgbm_oof", "morgan_plain", "prediction"),
    ("lgbm::mordred_plain", "lgbm_oof", "mordred_plain", "prediction"),
    ("lgbm::fp_mordred_plain", "lgbm_oof", "fp_mordred_plain", "prediction"),
    ("lgbm::external_class_plain", "lgbm_oof", "external_class_plain", "prediction"),
    ("lgbm::external_mixed_plain", "lgbm_oof", "external_mixed_plain", "prediction"),
    (
        "lgbm::external_class_mordred_plain",
        "lgbm_oof",
        "external_class_mordred_plain",
        "prediction",
    ),
    ("pfn::external_class_pca256", "pfn_oof", "external_class_pca256", "probability"),
    ("pfn::external_mixed_pca256", "pfn_oof", "external_mixed_pca256", "probability"),
    ("pfn::chemeleon_pca256", "pfn_oof", "chemeleon_pca256", "probability"),
    ("graph::tdi_only", "graph_oof", "tdi_only", "probability"),
    ("graph::tdi_aux", "graph_oof", "tdi_aux", "probability"),
]
atoms = {name: aligned(source, method, column) for name, source, method, column in specs}
if not np.allclose(np.mean(list(atoms.values()), axis=0), reference["prediction"].to_numpy()):
    raise ValueError("Atomic predictions do not reconstruct the reference old11 mean")
for extra in cfg.get("extra_candidates", []):
    atoms[extra["name"]] = aligned(extra["source"], extra["method"], extra["column"])

fixed_table = read("fixed_oof").filter(pl.col("endpoint") == endpoint)
fixed = {}
groups = None
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


def fit_caruana(x: np.ndarray, labels: np.ndarray, local_folds: np.ndarray, objective: str):
    counts = np.zeros(x.shape[1], dtype=int)
    total = np.zeros(len(labels), dtype=float)
    best_value, best_counts = -np.inf, None
    for step in range(1, cfg["caruana_steps"] + 1):
        values = [
            score_objective(labels, (total + x[:, index]) / step, local_folds, objective)
            for index in range(x.shape[1])
        ]
        choice = int(np.argmax(values))
        total += x[:, choice]
        counts[choice] += 1
        if values[choice] > best_value + 1e-12:
            best_value, best_counts = values[choice], counts.copy()
    return best_counts / best_counts.sum()


def crossfit_caruana(x: np.ndarray, objective: str):
    prediction = np.empty(len(y), dtype=float)
    weights = []
    for held in np.unique(folds):
        train, valid = folds != held, folds == held
        weight = fit_caruana(x[train], y[train], folds[train], objective)
        prediction[valid] = x[valid] @ weight
        weights.append({"held_fold": int(held), "weights": weight.tolist()})
    return prediction, weights


def crossfit_recipe(recipes: dict[str, np.ndarray], objective: str):
    prediction = np.empty(len(y), dtype=float)
    choices = []
    for held in np.unique(folds):
        train, valid = folds != held, folds == held
        values = {
            name: score_objective(y[train], pred[train], folds[train], objective)
            for name, pred in recipes.items()
        }
        choice = max(values, key=values.get)
        prediction[valid] = recipes[choice][valid]
        choices.append({"held_fold": int(held), "choice": choice, "value": values[choice]})
    return prediction, choices


predictions = {}
complexities = {}
selection_rows = []
member_names = list(atoms)
for variant, removed in cfg["ablations"].items():
    retained = [name for name in member_names if name not in removed]
    core = fold_rank(np.mean([atoms[name] for name in retained], axis=0), folds)
    live = (2 * core + fixed["logistic"]) / 3
    recipes = {
        "core": core,
        "core_logistic": live,
        "tail_rbf": (2 * live + fixed["rbf"]) / 3,
        "tail_complement": (2 * live + fixed["complement"]) / 3,
    }
    recipes["tail_average"] = (recipes["tail_rbf"] + recipes["tail_complement"]) / 2
    if variant == "full11":
        current = np.empty(len(y), dtype=float)
        for fold, recipe in enumerate(
            ["tail_rbf", "tail_rbf", "tail_complement", "tail_complement", "tail_complement"]
        ):
            mask = folds == fold
            current[mask] = recipes[recipe][mask]
        predictions["current_hierarchy"] = current
        complexities["current_hierarchy"] = 14
    for recipe, prediction in recipes.items():
        name = f"{variant}::{recipe}"
        predictions[name] = prediction
        complexities[name] = len(retained) + {
            "core": 0,
            "core_logistic": 1,
            "tail_rbf": 2,
            "tail_complement": 2,
            "tail_average": 3,
        }[recipe]
    for objective in ("auc", "mcc"):
        prediction, choices = crossfit_recipe(recipes, objective)
        name = f"{variant}::recipe_select_{objective}"
        predictions[name] = prediction
        complexities[name] = len(retained) + 3
        selection_rows.append({"model": name, "selection": json.dumps(choices)})
        blocks = np.column_stack([core, fixed["logistic"], fixed["rbf"], fixed["complement"]])
        prediction, weights = crossfit_caruana(blocks, objective)
        name = f"{variant}::block_caruana_{objective}"
        predictions[name] = prediction
        complexities[name] = len(retained) + 3
        selection_rows.append({"model": name, "selection": json.dumps(weights)})

# Search every non-empty subset of all 11 core members. Candidates are grouped
# by their actual number of component models, including downstream blocks. For
# each held fold the subset and recipe are chosen only on the other four folds.
candidate_groups = {}
recipe_cost = {
    "core": 0,
    "core_logistic": 1,
    "tail_rbf": 2,
    "tail_complement": 2,
    "tail_average": 3,
}
for size in range(1, len(member_names) + 1):
    for subset in combinations(member_names, size):
        core = fold_rank(np.mean([atoms[name] for name in subset], axis=0), folds)
        live = (2 * core + fixed["logistic"]) / 3
        recipes = {
            "core": core,
            "core_logistic": live,
            "tail_rbf": (2 * live + fixed["rbf"]) / 3,
            "tail_complement": (2 * live + fixed["complement"]) / 3,
        }
        recipes["tail_average"] = (recipes["tail_rbf"] + recipes["tail_complement"]) / 2
        for recipe, prediction in recipes.items():
            total = size + recipe_cost[recipe]
            candidate_groups.setdefault(total, []).append((subset, recipe, prediction))

for objective in ("auc", "mcc"):
    for total, candidates in sorted(candidate_groups.items()):
        prediction = np.empty(len(y), dtype=float)
        choices = []
        for held in np.unique(folds):
            train, valid = folds != held, folds == held
            values = [
                score_objective(y[train], candidate[2][train], folds[train], objective)
                for candidate in candidates
            ]
            choice = int(np.argmax(values))
            subset, recipe, candidate_prediction = candidates[choice]
            prediction[valid] = candidate_prediction[valid]
            choices.append(
                {
                    "held_fold": int(held),
                    "members": list(subset),
                    "recipe": recipe,
                    "training_objective": values[choice],
                }
            )
        full_values = [
            score_objective(y, candidate[2], folds, objective) for candidate in candidates
        ]
        full_choice = int(np.argmax(full_values))
        full_subset, full_recipe, _ = candidates[full_choice]
        choices.append(
            {
                "held_fold": "full",
                "members": list(full_subset),
                "recipe": full_recipe,
                "training_objective": full_values[full_choice],
            }
        )
        name = f"exhaustive_{objective}::components_{total}"
        predictions[name] = prediction
        complexities[name] = total
        selection_rows.append({"model": name, "selection": json.dumps(choices)})

rows, calls = [], {}
for name, prediction in predictions.items():
    selected = fixed_calls(prediction, folds)
    calls[name] = selected
    rows.append(
        {
            "model": name,
            "components": complexities[name],
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, selected)),
            "MCC_target_prior": expected_mcc(y, selected),
        }
    )

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
boot = []
for name in predictions:
    if name == "current_hierarchy":
        continue
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        sample_rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[sample_rows], calls[name][sample_rows])
            - expected_mcc(y[sample_rows], calls["current_hierarchy"][sample_rows])
        )
    values = np.asarray(deltas)
    boot.append(
        {
            "model": name,
            "mean_MCC_delta_vs_current": float(values.mean()),
            "ci025": float(np.quantile(values, 0.025)),
            "ci975": float(np.quantile(values, 0.975)),
            "probability_noninferior": float((values >= 0).mean()),
        }
    )

metrics = pl.DataFrame(rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(boot).write_csv(out / "bootstrap.csv")
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
(out / "summary.json").write_text(json.dumps({"selection": cfg["selection"]}, indent=2) + "\n")
print(metrics.head(25), flush=True)
