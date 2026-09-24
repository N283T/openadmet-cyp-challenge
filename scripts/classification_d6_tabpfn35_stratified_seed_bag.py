"""Evaluate label-blind support-selected stratified-random seed bags."""

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


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    calls[np.argsort(-values, kind="stable")[:count]] = True
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


base = (
    read("deployed_oof")
    .filter(pl.col("model") == "deployed_mean_block")
    .sort("compound_id")
)
ids = base.select("compound_id")
y = base["y_true"].to_numpy().astype(int)
deployed_oof = unit_rank(base["prediction"].to_numpy())
deployed_test_table = read("deployed_test").sort("compound_id")
test_ids = deployed_test_table.select("compound_id")
deployed_test = unit_rank(deployed_test_table["d6_score"].to_numpy())
oof_sources = pl.concat([read("primary_oof"), read("followup_oof")])
test_sources = pl.concat([read("primary_test"), read("followup_test")])


def aligned(
    table: pl.DataFrame, keys: pl.DataFrame, strategy: str, method: str
) -> tuple[np.ndarray, np.ndarray | None]:
    selected = table.filter(
        (pl.col("split_id") == strategy) & (pl.col("method") == method)
    ).select("compound_id", "prediction", *(["fold"] if "fold" in table.columns else []))
    joined = keys.join(selected, on="compound_id", validate="1:1", maintain_order="left")
    folds = joined["fold"].to_numpy().astype(int) if "fold" in joined.columns else None
    return joined["prediction"].to_numpy(), folds


predictions = {"deployed_mean_block": (deployed_oof, deployed_test)}
for method in cfg["methods"]:
    oof_members = []
    test_members = []
    for seed in cfg["selected_seeds"]:
        strategy = f"random_seed_{seed}"
        raw_oof, folds = aligned(oof_sources, ids, strategy, method)
        raw_test, _ = aligned(test_sources, test_ids, strategy, method)
        ranked_oof = fold_rank(raw_oof, folds)
        ranked_test = unit_rank(raw_test)
        predictions[f"{method}::seed_{seed}"] = (ranked_oof, ranked_test)
        oof_members.append(ranked_oof)
        test_members.append(ranked_test)
    predictions[f"{method}::support_selected_seed5_equal"] = (
        np.mean(oof_members, axis=0),
        np.mean(test_members, axis=0),
    )

base_test_calls = top_calls(deployed_test, cfg["test_positive_calls"])
metric_rows = []
oof_rows = []
test_rows = []
for name, (oof_prediction, test_prediction) in predictions.items():
    calls = top_calls(oof_prediction, cfg["oof_positive_calls"])
    test_calls = top_calls(test_prediction, cfg["test_positive_calls"])
    metric_rows.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, oof_prediction)),
            "AP": float(average_precision_score(y, oof_prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "test_swaps_from_deployed": int((test_calls & ~base_test_calls).sum()),
        }
    )
    oof_rows.append(
        ids.with_columns(
            pl.Series("y_true", y),
            pl.lit(name).alias("model"),
            pl.Series("prediction", oof_prediction),
            pl.Series("call", calls),
        )
    )
    test_rows.append(
        test_ids.with_columns(
            pl.lit(name).alias("model"),
            pl.Series("prediction", test_prediction),
            pl.Series("call", test_calls),
        )
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.concat(oof_rows).write_parquet(out / "oof.parquet")
pl.concat(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(json.dumps({"protocol": cfg["protocol"]}, indent=2) + "\n")
print(metrics, flush=True)
