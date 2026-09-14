"""Prespecified late fusion of endpoint-specific radius and Chemprop blocks."""

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


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        order = np.argsort(-scores[rows], kind="stable")
        calls[rows[order[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


compact_oof = read("compact_oof")
compact_test = read("compact_test")
radius_oof = read("radius_oof")
radius_test = read("radius_test")
metric_rows = []
oof_rows = []
test_rows = []
for endpoint in cfg["endpoints"]:
    settings = cfg["endpoint_settings"][endpoint]
    for strategy in cfg["strategies"]:
        radius = radius_oof.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == strategy)
            & (pl.col("method") == settings["radius_method"])
        ).sort("compound_id")
        chemprop = compact_oof.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == strategy)
            & (pl.col("method") == settings["chemprop_method"])
        ).sort("compound_id")
        if chemprop["compound_id"].to_list() != radius["compound_id"].to_list():
            raise ValueError(f"OOF keys differ for {endpoint} {strategy}")
        y = radius["y_true"].to_numpy().astype(int)
        folds = radius["fold"].to_numpy().astype(int)
        radius_rank = fold_rank(radius["prediction"].to_numpy(), folds)
        chemprop_rank = fold_rank(chemprop["prediction"].to_numpy(), folds)

        radius_deploy = radius_test.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == strategy)
            & (pl.col("method") == settings["radius_method"])
        ).sort("compound_id")
        chemprop_deploy = compact_test.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == strategy)
            & (pl.col("method") == settings["chemprop_method"])
        ).sort("compound_id")
        if chemprop_deploy["compound_id"].to_list() != radius_deploy["compound_id"].to_list():
            raise ValueError(f"Test keys differ for {endpoint} {strategy}")
        radius_test_rank = unit_rank(radius_deploy["prediction"].to_numpy())
        chemprop_test_rank = unit_rank(chemprop_deploy["prediction"].to_numpy())

        for method, weights in cfg["weights"].items():
            prediction = (weights[0] * radius_rank + weights[1] * chemprop_rank) / sum(weights)
            deploy_prediction = (
                weights[0] * radius_test_rank + weights[1] * chemprop_test_rank
            ) / sum(weights)
            calls = fixed_rate_calls(prediction, folds, settings["call_rate"])
            metric_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "method": method,
                    "radius_weight": weights[0],
                    "chemprop_weight": weights[1],
                    "AUROC": float(roc_auc_score(y, prediction)),
                    "AP": float(average_precision_score(y, prediction)),
                    "MCC_observed": float(matthews_corrcoef(y, calls)),
                    "MCC_target_prior": expected_mcc(y, calls, settings["target_prior"]),
                }
            )
            for row, value in enumerate(prediction):
                oof_rows.append(
                    {
                        "compound_id": radius["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": int(folds[row]),
                        "group": radius["group"][row],
                        "y_true": int(y[row]),
                        "method": method,
                        "prediction": float(value),
                        "call": bool(calls[row]),
                    }
                )
            for row, value in enumerate(deploy_prediction):
                test_rows.append(
                    {
                        "compound_id": radius_deploy["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "method": method,
                        "prediction": float(value),
                        "fold": -1,
                    }
                )

metrics = pl.DataFrame(metric_rows).sort(
    "endpoint", "strategy", "MCC_target_prior", descending=[False, False, True]
)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "scope": "Prespecified low-variance block fusion on both validation strategies.",
            "weights": cfg["weights"],
            "selection": "No candidate packaging or submission.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
