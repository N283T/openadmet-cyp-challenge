"""Distill the deployed six-member A4 Caruana ensemble into smaller subsets."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import minimize
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["member_oof"]).sort("compound_id")
test = pl.read_parquet(inp / cfg["member_test"]).sort("compound_id")
members = cfg["members"]
y = oof["y_true"].to_numpy().astype(int)
folds = oof["fold"].to_numpy().astype(int)
target_oof = oof[cfg["target"]].to_numpy()
target_test = test[cfg["target"]].to_numpy()
x_oof = np.column_stack([oof[name].to_numpy() for name in members])
x_test = np.column_stack([test[name].to_numpy() for name in members])


def fixed_rate_calls(values: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    result[np.argsort(-values, kind="stable")[:count]] = True
    return result


def fit_simplex(features: np.ndarray, target: np.ndarray) -> np.ndarray:
    if features.shape[1] == 1:
        return np.ones(1)
    initial = np.full(features.shape[1], 1 / features.shape[1])
    result = minimize(
        lambda weights: float(np.mean((features @ weights - target) ** 2)),
        initial,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * features.shape[1],
        constraints={"type": "eq", "fun": lambda weights: weights.sum() - 1.0},
        options={"ftol": 1e-14, "maxiter": 2000},
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = np.clip(result.x, 0, 1)
    return weights / weights.sum()


target_oof_calls = fixed_rate_calls(target_oof)
target_test_calls = top_calls(target_test, cfg["test_positive_calls"])
rows = []
predictions = {}
for subset_size in range(1, len(members) + 1):
    for subset in itertools.combinations(range(len(members)), subset_size):
        indices = list(subset)
        weights = fit_simplex(x_oof[:, indices], target_oof)
        prediction = x_oof[:, indices] @ weights
        deploy_prediction = x_test[:, indices] @ weights
        oof_calls = fixed_rate_calls(prediction)
        test_calls = top_calls(deploy_prediction, cfg["test_positive_calls"])
        names = [members[index] for index in subset]
        key = " + ".join(names)
        predictions[key] = deploy_prediction
        rows.append(
            {
                "requested_size": subset_size,
                "effective_size": int((weights > cfg["effective_weight_min"]).sum()),
                "members": key,
                "weights": json.dumps(dict(zip(names, weights.tolist(), strict=True))),
                "oof_rmse_to_original": float(np.sqrt(np.mean((prediction - target_oof) ** 2))),
                "oof_spearman_to_original": float(spearmanr(prediction, target_oof).statistic),
                "oof_auroc": float(roc_auc_score(y, prediction)),
                "oof_mcc_at_fold_rate": float(matthews_corrcoef(y, oof_calls)),
                "oof_call_swaps_each_way": int(np.logical_and(oof_calls, ~target_oof_calls).sum()),
                "test_rmse_to_original": float(
                    np.sqrt(np.mean((deploy_prediction - target_test) ** 2))
                ),
                "test_spearman_to_original": float(
                    spearmanr(deploy_prediction, target_test).statistic
                ),
                "test_call_swaps_each_way": int(
                    np.logical_and(test_calls, ~target_test_calls).sum()
                ),
            }
        )

results = pl.DataFrame(rows)
results.write_csv(out / "all-subsets.csv")
best = (
    results.sort("requested_size", "oof_rmse_to_original")
    .group_by("requested_size", maintain_order=True)
    .first()
    .sort("requested_size")
)
best.write_csv(out / "best-by-size.csv")
prediction_table = test.select("compound_id", cfg["target"])
for row in best.iter_rows(named=True):
    prediction_table = prediction_table.with_columns(
        pl.Series(f"compressed_size_{row['requested_size']}", predictions[row["members"]])
    )
prediction_table.write_parquet(out / "best-test-predictions.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "objective": "OOF probability RMSE to the deployed full-OOF Caruana mixture",
            "selection": "Best subset at each size selected without blind labels",
            "test_positive_calls": cfg["test_positive_calls"],
            "caveat": (
                "This is prediction distillation, not an independently selected performance "
                "improvement. OOF metrics are descriptive because the original deployment "
                "weights were selected on full OOF labels."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(best, flush=True)
