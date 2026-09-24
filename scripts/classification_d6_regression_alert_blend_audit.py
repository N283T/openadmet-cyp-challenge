"""Audit fixed rank additions to the currently deployed CYP2D6 score."""

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
endpoint = cfg["endpoint"]


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    ranked = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        ranked[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return ranked


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


# Reconstruct the cross-fitted mean-block OOF prediction used to justify the
# deployed full-data test recipe. This is the exact logic of run b3a7120e.
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
    if not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete alignment for {source}/{method}")
    return joined["prediction"].to_numpy()


specs = {
    "lgbm::morgan_plain": ("lgbm_oof", "morgan_plain", "prediction"),
    "lgbm::mordred_plain": ("lgbm_oof", "mordred_plain", "prediction"),
    "lgbm::fp_mordred_plain": ("lgbm_oof", "fp_mordred_plain", "prediction"),
    "lgbm::external_class_plain": ("lgbm_oof", "external_class_plain", "prediction"),
    "lgbm::external_mixed_plain": ("lgbm_oof", "external_mixed_plain", "prediction"),
    "lgbm::external_class_mordred_plain": (
        "lgbm_oof",
        "external_class_mordred_plain",
        "prediction",
    ),
    "pfn::external_class_pca256": ("pfn_oof", "external_class_pca256", "probability"),
    "pfn::external_mixed_pca256": ("pfn_oof", "external_mixed_pca256", "probability"),
    "pfn::chemeleon_pca256": ("pfn_oof", "chemeleon_pca256", "probability"),
    "graph::tdi_only": ("graph_oof", "tdi_only", "probability"),
    "graph::tdi_aux": ("graph_oof", "tdi_aux", "probability"),
    "pfn::external_mixed_raw600": (
        "raw600_oof",
        "external_mixed_raw600",
        "probability",
    ),
}
atoms = {name: aligned(*source) for name, source in specs.items()}
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


def recipe_prediction(members: list[str], recipe: str) -> np.ndarray:
    values = []
    for member in members:
        if member == cfg["pca_member"]:
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
choices = [
    choice
    for choice in json.loads(selection["selection"])
    if choice["held_fold"] != "full"
]
base_oof = np.empty(len(y), dtype=float)
for choice in choices:
    candidate = recipe_prediction(choice["members"], choice["recipe"])
    mask = folds == int(choice["held_fold"])
    base_oof[mask] = candidate[mask]
base_oof = fold_rank(base_oof, folds)
base_calls = fixed_calls(base_oof, folds)

deployed_test = read("deployed_test").sort("compound_id")
test_ids = deployed_test.select("compound_id")
base_test = unit_rank(deployed_test["d6_score"].to_numpy())
base_test_calls = top_calls(base_test)
if not np.array_equal(base_test_calls, deployed_test["d6_call"].to_numpy().astype(bool)):
    raise ValueError("Deployed D6 test score does not reconstruct its recorded calls")


def candidate_values(spec: dict, partition: str) -> np.ndarray:
    source = f"{spec['source']}_{partition}"
    table = (
        read(source)
        .filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("split_id") == cfg["split_id"])
            & (pl.col("method") == spec["method"])
        )
        .select("compound_id", "prediction")
    )
    keys = ids if partition == "oof" else test_ids
    joined = keys.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete candidate {spec['name']} {partition}")
    values = joined["prediction"].to_numpy()
    return fold_rank(values, folds) if partition == "oof" else unit_rank(values)


metric_rows = []
test_rows = []
bootstrap_rows = []
oof_rows = []


def add_metric(name: str, prediction: np.ndarray, calls: np.ndarray) -> None:
    metric_rows.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "positive_calls": int(calls.sum()),
            "swaps_from_base": int((calls & ~base_calls).sum()),
        }
    )
    for row in range(len(y)):
        oof_rows.append(
            {
                "compound_id": ids["compound_id"][row],
                "fold": int(folds[row]),
                "group": groups[row],
                "y_true": int(y[row]),
                "model": name,
                "prediction": float(prediction[row]),
                "call": bool(calls[row]),
            }
        )


add_metric("deployed_mean_block", base_oof, base_calls)
test_rows.append(
    {
        "model": "deployed_mean_block",
        "positive_calls": int(base_test_calls.sum()),
        "swaps_from_deployed": 0,
    }
)

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
for spec in cfg["candidates"]:
    candidate_oof = candidate_values(spec, "oof")
    candidate_test = candidate_values(spec, "test")
    for weight in cfg["candidate_weights"]:
        name = f"{spec['name']}::w{weight:g}"
        prediction = (1 - weight) * base_oof + weight * candidate_oof
        calls = fixed_calls(prediction, folds)
        add_metric(name, prediction, calls)
        test_prediction = (1 - weight) * base_test + weight * candidate_test
        test_calls = top_calls(test_prediction)
        test_rows.append(
            {
                "model": name,
                "positive_calls": int(test_calls.sum()),
                "swaps_from_deployed": int((test_calls & ~base_test_calls).sum()),
            }
        )
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows = np.concatenate([group_rows[group] for group in sampled])
            deltas.append(
                expected_mcc(y[rows], calls[rows])
                - expected_mcc(y[rows], base_calls[rows])
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
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": endpoint,
            "selection": cfg["selection"],
            "base_OOF_reconstruction": {
                "AUROC": float(roc_auc_score(y, base_oof)),
                "AP": float(average_precision_score(y, base_oof)),
                "MCC_target_prior": expected_mcc(y, base_calls),
            },
            "test_policy": "Keep 49 highest blended ranks; CYP3A4 is untouched.",
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
