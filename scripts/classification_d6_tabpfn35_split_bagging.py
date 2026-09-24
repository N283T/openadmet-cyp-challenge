"""Audit equal split-bagging for standalone D6 TabPFN 3.5 models."""

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


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    calls[np.argsort(-values, kind="stable")[:count]] = True
    return calls


def weighted_calls(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    order = np.argsort(-values, kind="stable")
    target = cfg["oof_positive_calls"] / len(values) * weights.sum()
    count = int(np.searchsorted(np.cumsum(weights[order]), target, side="left") + 1)
    calls = np.zeros(len(values), dtype=bool)
    calls[order[:count]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, weights: np.ndarray) -> float:
    positive = y == 1
    tpr = float(np.average(calls[positive], weights=weights[positive]))
    fpr = float(np.average(calls[~positive], weights=weights[~positive]))
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


base = (
    read("deployed_oof")
    .filter(pl.col("model") == cfg["deployed_model"])
    .sort("compound_id")
)
ids = base.select("compound_id")
y = base["y_true"].to_numpy().astype(int)
groups = base["group"].cast(pl.String).to_numpy()
deployed_oof = unit_rank(base["prediction"].to_numpy())
support = ids.join(
    read("support").select("compound_id", "importance_weight"),
    on="compound_id",
    validate="1:1",
    maintain_order="left",
)
weights = support["importance_weight"].to_numpy()


def aligned(source: str, split_id: str, method: str) -> tuple[np.ndarray, np.ndarray]:
    table = (
        read(source)
        .filter((pl.col("split_id") == split_id) & (pl.col("method") == method))
        .select("compound_id", "fold", "y_true", "prediction")
    )
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if not np.array_equal(joined["y_true"].to_numpy(), y):
        raise ValueError(f"Label mismatch for {source}/{split_id}/{method}")
    folds = joined["fold"].to_numpy().astype(int)
    return fold_rank(joined["prediction"].to_numpy(), folds), folds


def aligned_test(source: str, split_id: str, method: str, test_ids: pl.DataFrame) -> np.ndarray:
    table = (
        read(source)
        .filter((pl.col("split_id") == split_id) & (pl.col("method") == method))
        .select("compound_id", "prediction")
    )
    joined = test_ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    return unit_rank(joined["prediction"].to_numpy())


deployed_test_table = read("deployed_test").sort("compound_id")
test_ids = deployed_test_table.select("compound_id")
deployed_test = unit_rank(deployed_test_table["d6_score"].to_numpy())

predictions: dict[str, tuple[np.ndarray, np.ndarray]] = {
    "deployed_mean_block": (deployed_oof, deployed_test)
}
for method in cfg["methods"]:
    umap, _ = aligned("dual_oof", "umap", method)
    scaffold, _ = aligned("dual_oof", "scaffold", method)
    random, _ = aligned("random_oof", "random", method)
    umap_test = aligned_test("dual_test", "umap", method, test_ids)
    scaffold_test = aligned_test("dual_test", "scaffold", method, test_ids)
    random_test = aligned_test("random_test", "random", method, test_ids)
    prefix = f"tabpfn35fast::{method}"
    predictions[f"{prefix}::umap"] = (umap, umap_test)
    predictions[f"{prefix}::scaffold"] = (scaffold, scaffold_test)
    predictions[f"{prefix}::random"] = (random, random_test)
    predictions[f"{prefix}::umap_random_equal"] = (
        (umap + random) / 2,
        (umap_test + random_test) / 2,
    )
    predictions[f"{prefix}::all3_equal"] = (
        (umap + scaffold + random) / 3,
        (umap_test + scaffold_test + random_test) / 3,
    )

metric_rows = []
oof_rows = []
test_rows = []
calls_by_model = {}
for name, (oof_prediction, test_prediction) in predictions.items():
    calls = top_calls(oof_prediction, cfg["oof_positive_calls"])
    weighted = weighted_calls(oof_prediction, weights)
    test_calls = top_calls(test_prediction, cfg["test_positive_calls"])
    calls_by_model[name] = calls
    metric_rows.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, oof_prediction)),
            "AP": float(average_precision_score(y, oof_prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls, np.ones(len(y))),
            "support_weighted_AUROC": float(
                roc_auc_score(y, oof_prediction, sample_weight=weights)
            ),
            "support_weighted_AP": float(
                average_precision_score(y, oof_prediction, sample_weight=weights)
            ),
            "support_weighted_MCC": expected_mcc(y, weighted, weights),
            "test_swaps_from_deployed": int((test_calls & ~top_calls(deployed_test, 49)).sum()),
        }
    )
    oof_rows.append(
        ids.with_columns(
            pl.Series("y_true", y),
            pl.lit(name).alias("model"),
            pl.Series("prediction", oof_prediction),
            pl.Series("call", calls),
            pl.Series("support_weighted_call", weighted),
        )
    )
    test_rows.append(
        test_ids.with_columns(
            pl.lit(name).alias("model"),
            pl.Series("prediction", test_prediction),
            pl.Series("call", test_calls),
        )
    )

rng = np.random.default_rng(cfg["bootstrap_seed"])
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
bootstrap_rows = []
for name, calls in calls_by_model.items():
    if name == "deployed_mean_block":
        continue
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], calls[rows], np.ones(len(rows)))
            - expected_mcc(
                y[rows], calls_by_model["deployed_mean_block"][rows], np.ones(len(rows))
            )
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
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.concat(oof_rows).write_parquet(out / "oof.parquet")
pl.concat(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(json.dumps({"protocol": cfg["protocol"]}, indent=2) + "\n")
print(metrics, flush=True)
