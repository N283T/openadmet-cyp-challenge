"""Equal-vote D6 consensus across complete TabPFN split grids."""

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


oof_long = read("split_bag_oof")
test_long = read("split_bag_test")
base = oof_long.filter(pl.col("model") == "deployed_mean_block").sort("compound_id")
ids = base.select("compound_id")
y = base["y_true"].to_numpy().astype(int)
support = ids.join(
    read("support").select("compound_id", "importance_weight"),
    on="compound_id",
    validate="1:1",
    maintain_order="left",
)
weights = support["importance_weight"].to_numpy()
test_base = test_long.filter(pl.col("model") == "deployed_mean_block").sort("compound_id")
test_ids = test_base.select("compound_id")


def aligned(table: pl.DataFrame, keys: pl.DataFrame, model: str) -> np.ndarray:
    selected = table.filter(pl.col("model") == model).select("compound_id", "prediction")
    return keys.join(selected, on="compound_id", validate="1:1", maintain_order="left")[
        "prediction"
    ].to_numpy()


predictions: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
base_oof = aligned(oof_long, ids, "deployed_mean_block")
base_test = aligned(test_long, test_ids, "deployed_mean_block")
predictions["deployed_mean_block"] = (
    base_oof,
    base_test,
    top_calls(base_oof, cfg["oof_positive_calls"]),
    top_calls(base_test, cfg["test_positive_calls"]),
)

for grid_name, suffixes in [
    ("support6", cfg["support6_suffixes"]),
    ("all9", cfg["all9_suffixes"]),
]:
    members = [f"{root}::{suffix}" for root in cfg["method_roots"] for suffix in suffixes]
    oof_matrix = np.column_stack([unit_rank(aligned(oof_long, ids, name)) for name in members])
    test_matrix = np.column_stack(
        [unit_rank(aligned(test_long, test_ids, name)) for name in members]
    )
    oof_mean = oof_matrix.mean(axis=1)
    test_mean = test_matrix.mean(axis=1)
    oof_member_calls = np.column_stack(
        [
            top_calls(oof_matrix[:, index], cfg["oof_positive_calls"])
            for index in range(len(members))
        ]
    )
    test_member_calls = np.column_stack(
        [
            top_calls(test_matrix[:, index], cfg["test_positive_calls"])
            for index in range(len(members))
        ]
    )
    oof_vote_score = oof_member_calls.sum(axis=1) + oof_mean
    test_vote_score = test_member_calls.sum(axis=1) + test_mean
    predictions[f"{grid_name}_equal_rank"] = (
        oof_mean,
        test_mean,
        top_calls(oof_mean, cfg["oof_positive_calls"]),
        top_calls(test_mean, cfg["test_positive_calls"]),
    )
    predictions[f"{grid_name}_top_call_vote"] = (
        oof_vote_score,
        test_vote_score,
        top_calls(oof_vote_score, cfg["oof_positive_calls"]),
        top_calls(test_vote_score, cfg["test_positive_calls"]),
    )

base_test_calls = predictions["deployed_mean_block"][3]
metric_rows = []
oof_rows = []
test_rows = []
for name, (oof_score, test_score, calls, test_calls) in predictions.items():
    weighted = weighted_calls(oof_score, weights)
    metric_rows.append(
        {
            "model": name,
            "AUROC": float(roc_auc_score(y, oof_score)),
            "AP": float(average_precision_score(y, oof_score)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls, np.ones(len(y))),
            "support_weighted_AUROC": float(roc_auc_score(y, oof_score, sample_weight=weights)),
            "support_weighted_AP": float(
                average_precision_score(y, oof_score, sample_weight=weights)
            ),
            "support_weighted_MCC": expected_mcc(y, weighted, weights),
            "test_swaps_from_deployed": int((test_calls & ~base_test_calls).sum()),
        }
    )
    oof_rows.append(
        ids.with_columns(
            pl.Series("y_true", y),
            pl.lit(name).alias("model"),
            pl.Series("score", oof_score),
            pl.Series("call", calls),
            pl.Series("support_weighted_call", weighted),
        )
    )
    test_rows.append(
        test_ids.with_columns(
            pl.lit(name).alias("model"),
            pl.Series("score", test_score),
            pl.Series("call", test_calls),
        )
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.concat(oof_rows).write_parquet(out / "oof.parquet")
pl.concat(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(json.dumps({"protocol": cfg["protocol"]}, indent=2) + "\n")
print(metrics, flush=True)
