"""Compare A4 raw replacements and simplify a common prediction bank."""

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
endpoint = cfg["endpoint"]


def read(key):
    return pl.read_parquet(inp / cfg[key])


def unit_rank(values):
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values, folds):
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_calls(values, folds):
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def top_calls(values):
    result = np.zeros(len(values), dtype=bool)
    result[np.argsort(-values, kind="stable")[: cfg["test_calls"]]] = True
    return result


def expected_mcc(y, calls):
    positive = y == 1
    tpr, fpr = float(calls[positive].mean()), float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def objective(y, prediction, folds, name):
    if name == "auc":
        return float(roc_auc_score(y, prediction))
    return expected_mcc(y, fixed_calls(prediction, folds))


def align_long(table, ids, method, test=False):
    predicate = (pl.col("endpoint") == endpoint) & (pl.col("method") == method)
    if "strategy" in table.columns:
        predicate &= pl.col("strategy") == "umap"
    column = "probability" if "probability" in table.columns else "prediction"
    selected = table.filter(predicate)
    if test and selected.select("compound_id").is_duplicated().any():
        selected = selected.group_by("compound_id").agg(pl.col(column).mean())
    selected = selected.select("compound_id", pl.col(column).alias("prediction"))
    joined = ids.join(selected, on="compound_id", validate="1:1", maintain_order="left")
    if joined.height != ids.height or not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete alignment for {method}")
    return joined["prediction"].to_numpy()


member_oof = read("member_oof").sort("compound_id")
member_test = read("member_test").sort("compound_id")
ids, test_ids = member_oof.select("compound_id"), member_test.select("compound_id")
folds = member_oof["fold"].to_numpy().astype(int)
y = member_oof["y_true"].to_numpy().astype(int)

names = [
    "historical",
    "jeremy_tabicl_pca256",
    "jeremy_tabicl_raw600",
    "jeremy_tabpfn_pca256",
    "jeremy_tabpfn_raw600",
    "external_class_mordred",
    "pxr_admetica",
    "external_class_raw600",
    "external_mixed_raw600",
    "chemeleon_parent_raw_lgbm",
]
member_columns = {
    "historical": "history::historical_stack_mordred",
    "jeremy_tabicl_pca256": "jeremy8::tabicl_jeremy_pca256",
    "jeremy_tabpfn_pca256": "jeremy8::tabpfn_jeremy_pca256",
    "external_class_mordred": "lgbm::external_class_mordred_balanced",
}
oof, test = {}, {}
for name, column in member_columns.items():
    oof[name] = fold_rank(member_oof[column].to_numpy(), folds)
    test[name] = unit_rank(member_test[column].to_numpy())

raw_oof, raw_test = read("raw_oof"), read("raw_test")
for name, method in {
    "jeremy_tabicl_raw600": "tabicl_jeremy_raw600",
    "jeremy_tabpfn_raw600": "tabpfn_jeremy_raw600",
    "external_class_raw600": "tabpfn_external_class_raw600",
    "external_mixed_raw600": "tabpfn_external_mixed_raw600",
}.items():
    oof[name] = fold_rank(align_long(raw_oof, ids, method), folds)
    test[name] = unit_rank(align_long(raw_test, test_ids, method, True))

oof["pxr_admetica"] = fold_rank(align_long(read("aux_oof"), ids, cfg["aux_method"]), folds)
test["pxr_admetica"] = unit_rank(align_long(read("aux_test"), test_ids, cfg["aux_method"], True))
oof["chemeleon_parent_raw_lgbm"] = fold_rank(
    align_long(read("embedding_oof"), ids, cfg["embedding_method"]), folds
)
test["chemeleon_parent_raw_lgbm"] = unit_rank(
    align_long(read("embedding_test"), test_ids, cfg["embedding_method"], True)
)

current = 0.8 * fold_rank(member_oof["caruana_auc"].to_numpy(), folds) + 0.2 * oof["pxr_admetica"]
current_test = 0.8 * unit_rank(member_test["caruana_auc"].to_numpy()) + 0.2 * test["pxr_admetica"]

member_metrics = []
for name in names:
    selected = fixed_calls(oof[name], folds)
    member_metrics.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, oof[name])),
            "AP": float(average_precision_score(y, oof[name])),
            "MCC_observed": float(matthews_corrcoef(y, selected)),
            "MCC_target_prior": expected_mcc(y, selected),
        }
    )


def fit_caruana(matrix, labels, local_folds, objective_name):
    counts = np.zeros(matrix.shape[1], dtype=int)
    total = np.zeros(len(labels), dtype=float)
    best_value, best_counts = -np.inf, None
    for step in range(1, cfg["caruana_steps"] + 1):
        values = [
            objective(
                labels,
                (total + matrix[:, index]) / step,
                local_folds,
                objective_name,
            )
            for index in range(matrix.shape[1])
        ]
        choice = int(np.argmax(values))
        total += matrix[:, choice]
        counts[choice] += 1
        if values[choice] > best_value + 1e-12:
            best_value, best_counts = values[choice], counts.copy()
    return best_counts / best_counts.sum()


predictions = {"current_fixed_80_20": current}
test_predictions = {"current_fixed_80_20": current_test}
selection_rows = []
for pool_name, pool_names in cfg["pools"].items():
    matrix = np.column_stack([oof[name] for name in pool_names])
    test_matrix = np.column_stack([test[name] for name in pool_names])
    for objective_name in ["auc", "mcc"]:
        prediction = np.empty(len(y), dtype=float)
        for held in np.unique(folds):
            train, valid = folds != held, folds == held
            weights = fit_caruana(matrix[train], y[train], folds[train], objective_name)
            prediction[valid] = matrix[valid] @ weights
            selection_rows.append(
                {
                    "model": f"{pool_name}::caruana_{objective_name}",
                    "held_fold": int(held),
                    "members": json.dumps(pool_names),
                    "weights": json.dumps(dict(zip(pool_names, weights.tolist(), strict=True))),
                }
            )
        weights = fit_caruana(matrix, y, folds, objective_name)
        model = f"{pool_name}::caruana_{objective_name}"
        predictions[model] = prediction
        test_predictions[model] = test_matrix @ weights
        selection_rows.append(
            {
                "model": model,
                "held_fold": -1,
                "members": json.dumps(pool_names),
                "weights": json.dumps(dict(zip(pool_names, weights.tolist(), strict=True))),
            }
        )

# Exhaustive equal-weight subsets keep component count and deployment transparent.
for objective_name in ["auc", "mcc"]:
    for size in range(1, len(names) + 1):
        subsets = list(combinations(names, size))
        subset_predictions = [np.mean([oof[name] for name in subset], axis=0) for subset in subsets]
        prediction = np.empty(len(y), dtype=float)
        choices = []
        for held in np.unique(folds):
            train, valid = folds != held, folds == held
            values = [
                objective(y[train], candidate[train], folds[train], objective_name)
                for candidate in subset_predictions
            ]
            choice = int(np.argmax(values))
            prediction[valid] = subset_predictions[choice][valid]
            choices.append({"held_fold": int(held), "members": list(subsets[choice])})
        full_values = [
            objective(y, candidate, folds, objective_name) for candidate in subset_predictions
        ]
        full_choice = int(np.argmax(full_values))
        full_subset = subsets[full_choice]
        model = f"exhaustive_{objective_name}::components_{size}"
        predictions[model] = prediction
        test_predictions[model] = np.mean([test[name] for name in full_subset], axis=0)
        choices.append({"held_fold": -1, "members": list(full_subset)})
        selection_rows.append(
            {
                "model": model,
                "held_fold": -2,
                "members": json.dumps(choices),
                "weights": None,
            }
        )

metric_rows, calls = [], {}
for model, prediction in predictions.items():
    selected = fixed_calls(prediction, folds)
    calls[model] = selected
    metric_rows.append(
        {
            "model": model,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, selected)),
            "MCC_target_prior": expected_mcc(y, selected),
            "test_call_swaps_vs_current": int(
                (top_calls(test_predictions[model]) & ~top_calls(current_test)).sum()
            ),
        }
    )

group_source = (
    read("embedding_oof")
    .filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("strategy") == "umap")
        & (pl.col("method") == cfg["embedding_method"])
    )
    .select("compound_id", pl.col("group").cast(pl.String))
)
groups = ids.join(group_source, on="compound_id", validate="1:1", maintain_order="left")[
    "group"
].to_numpy()
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
bootstrap_rows = []
for model in predictions:
    if model == "current_fixed_80_20":
        continue
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], calls[model][rows])
            - expected_mcc(y[rows], calls["current_fixed_80_20"][rows])
        )
    delta = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "model": model,
            "mean_MCC_delta": float(delta.mean()),
            "ci025": float(np.quantile(delta, 0.025)),
            "ci975": float(np.quantile(delta, 0.975)),
            "probability_positive": float((delta > 0).mean()),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(member_metrics).sort("MCC_target_prior", descending=True).write_csv(
    out / "member-metrics.csv"
)
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
pl.concat(
    [
        ids.with_columns(
            pl.lit(model).alias("model"),
            pl.Series("fold", folds),
            pl.Series("y_true", y),
            pl.Series("prediction", prediction, dtype=pl.Float64),
            pl.Series("call", calls[model]),
        )
        for model, prediction in predictions.items()
    ]
).write_parquet(out / "oof.parquet")
pl.concat(
    [
        test_ids.with_columns(
            pl.lit(model).alias("model"),
            pl.Series("prediction", prediction, dtype=pl.Float64),
            pl.Series("call", top_calls(prediction)),
        )
        for model, prediction in test_predictions.items()
    ]
).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(json.dumps({"selection": cfg["selection"]}, indent=2) + "\n")
print(metrics.head(30), flush=True)
