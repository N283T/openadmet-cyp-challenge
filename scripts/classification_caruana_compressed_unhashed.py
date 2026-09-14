"""Test an unhashed-radius block inside the compressed A4 Caruana recipe."""

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


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_calls(values: np.ndarray, folds: np.ndarray | None = None) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    if folds is None:
        result[np.argsort(-values, kind="stable")[: cfg["positive_calls"]]] = True
        return result
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        order = np.argsort(-values[positions], kind="stable")
        result[positions[order[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


members_oof = pl.read_parquet(inp / cfg["member_oof"]).sort("compound_id")
members_test = pl.read_parquet(inp / cfg["member_test"]).sort("compound_id")
unhashed_methods = (
    cfg["unhashed_methods"] if "unhashed_methods" in cfg else [cfg["unhashed_method"]]
)
radius_oof_all = pl.read_parquet(inp / cfg["radius_oof"]).filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("strategy") == cfg["strategy"])
)
radius_test_all = pl.read_parquet(inp / cfg["radius_test"]).filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("strategy") == cfg["strategy"])
)
radius_oof = {
    method: radius_oof_all.filter(pl.col("method") == method).sort("compound_id")
    for method in unhashed_methods
}
radius_test = {
    method: radius_test_all.filter(pl.col("method") == method).sort("compound_id")
    for method in unhashed_methods
}
keys = ["compound_id", "fold", "y_true"]
for method in unhashed_methods:
    if not members_oof.select(keys).equals(radius_oof[method].select(keys)):
        raise ValueError(f"Member and {method} OOF keys differ")
    if members_test["compound_id"].to_list() != radius_test[method]["compound_id"].to_list():
        raise ValueError(f"Member and {method} test IDs differ")

y = members_oof["y_true"].to_numpy().astype(int)
folds = members_oof["fold"].to_numpy().astype(int)
groups = radius_oof[unhashed_methods[0]]["group"].cast(pl.String).to_numpy()
compressed_oof = sum(
    weight * members_oof[member].to_numpy()
    for member, weight in zip(
        [cfg["historical_member"], cfg["jeremy_member"]],
        cfg["compressed_weights"],
        strict=True,
    )
)
compressed_test = sum(
    weight * members_test[member].to_numpy()
    for member, weight in zip(
        [cfg["historical_member"], cfg["jeremy_member"]],
        cfg["compressed_weights"],
        strict=True,
    )
)

oof_predictions = {
    "full_caruana": fold_rank(members_oof["caruana_auc"].to_numpy(), folds),
    "compressed2": fold_rank(compressed_oof, folds),
}
test_predictions = {
    "full_caruana": unit_rank(members_test["caruana_auc"].to_numpy()),
    "compressed2": unit_rank(compressed_test),
}
for method in unhashed_methods:
    short_name = (
        "pure_unhashed" if method == "radius_specific_unhashed_count" else "morgan_unhashed"
    )
    block_name = f"block_{short_name}"
    oof_predictions[block_name] = fold_rank(radius_oof[method]["prediction"].to_numpy(), folds)
    test_predictions[block_name] = unit_rank(radius_test[method]["prediction"].to_numpy())
    for base_weight, unhashed_weight in cfg["blend_weights"]:
        name = f"compressed2_{base_weight}_{short_name}_{unhashed_weight}"
        total = base_weight + unhashed_weight
        oof_predictions[name] = (
            base_weight * oof_predictions["compressed2"]
            + unhashed_weight * oof_predictions[block_name]
        ) / total
        test_predictions[name] = (
            base_weight * test_predictions["compressed2"]
            + unhashed_weight * test_predictions[block_name]
        ) / total

baseline_oof_calls = fixed_calls(oof_predictions["full_caruana"], folds)
compressed_oof_calls = fixed_calls(oof_predictions["compressed2"], folds)
baseline_test_calls = fixed_calls(test_predictions["full_caruana"])
compressed_test_calls = fixed_calls(test_predictions["compressed2"])
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
metric_rows = []
bootstrap_rows = []
test_rows = []
for name, prediction in oof_predictions.items():
    calls = fixed_calls(prediction, folds)
    test_prediction = test_predictions[name]
    test_calls = fixed_calls(test_prediction)
    metric_rows.append(
        {
            "method": name,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "test_spearman_to_full": float(
                spearmanr(test_prediction, test_predictions["full_caruana"]).statistic
            ),
            "test_swaps_vs_full_each_way": int((test_calls & ~baseline_test_calls).sum()),
            "test_swaps_vs_compressed2_each_way": int((test_calls & ~compressed_test_calls).sum()),
        }
    )
    if name.startswith("compressed2_") and name != "compressed2":
        deltas_vs_full = []
        deltas_vs_compressed = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows = np.concatenate([group_rows[group] for group in sampled])
            score = expected_mcc(y[rows], calls[rows])
            deltas_vs_full.append(score - expected_mcc(y[rows], baseline_oof_calls[rows]))
            deltas_vs_compressed.append(score - expected_mcc(y[rows], compressed_oof_calls[rows]))
        bootstrap_rows.append(
            {
                "method": name,
                "delta_vs_full": expected_mcc(y, calls) - expected_mcc(y, baseline_oof_calls),
                "delta_vs_compressed2": expected_mcc(y, calls)
                - expected_mcc(y, compressed_oof_calls),
                "delta_vs_full_p025": float(np.quantile(deltas_vs_full, 0.025)),
                "delta_vs_full_p975": float(np.quantile(deltas_vs_full, 0.975)),
                "delta_vs_compressed2_p025": float(np.quantile(deltas_vs_compressed, 0.025)),
                "delta_vs_compressed2_p975": float(np.quantile(deltas_vs_compressed, 0.975)),
            }
        )
    for row, value in enumerate(test_prediction):
        test_rows.append(
            {
                "compound_id": members_test["compound_id"][row],
                "endpoint": cfg["endpoint"],
                "method": name,
                "prediction": float(value),
                "call": bool(test_calls[row]),
                "fold": -1,
            }
        )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).sort("delta_vs_full", descending=True).write_csv(
    out / "conditional-bootstrap.csv"
)
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "scope": "Exploratory fixed-weight addition to compressed A4 Caruana",
            "unhashed_methods": unhashed_methods,
            "selection_warning": "The small blend grid is compared on the same OOF labels.",
            "submission": "None",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
print(pl.DataFrame(bootstrap_rows), flush=True)
