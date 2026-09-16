"""Paired group bootstrap for fixed CYP-MAP mechanism readout predictions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["oof"]).filter(pl.col("endpoint") == cfg["endpoint"])
rng = np.random.default_rng(cfg["seed"])


def projected_mcc(y, calls, weights, prevalence):
    positive = y == 1
    tpr = np.average(calls[positive], weights=weights[positive])
    fpr = np.average(calls[~positive], weights=weights[~positive])
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


rows = []
for strategy in cfg["strategies"]:
    frame = oof.filter(pl.col("strategy") == strategy)
    anchor = frame.filter(pl.col("method") == cfg["anchor"]).sort("compound_id")
    groups = anchor["group"].to_numpy()
    unique_groups = np.unique(groups)
    y = anchor["y_true"].to_numpy().astype(int)
    anchor_scores = anchor["prediction"].to_numpy()
    anchor_calls = anchor["call"].to_numpy()
    for method in cfg["methods"]:
        candidate = frame.filter(pl.col("method") == method).sort("compound_id")
        if not candidate["compound_id"].equals(anchor["compound_id"]):
            raise ValueError(f"OOF alignment failed for {method} {strategy}")
        scores = candidate["prediction"].to_numpy()
        calls = candidate["call"].to_numpy()
        draws = {"AUROC": [], "AP": [], "MCC_target_prior": []}
        for _ in range(cfg["replicates"]):
            sampled = rng.choice(unique_groups, size=len(unique_groups), replace=True)
            sampled_groups, sampled_counts = np.unique(sampled, return_counts=True)
            counts = dict(zip(sampled_groups, sampled_counts, strict=True))
            weights = np.asarray([counts.get(group, 0) for group in groups], dtype=float)
            if np.unique(y[weights > 0]).size < 2:
                continue
            draws["AUROC"].append(
                roc_auc_score(y, scores, sample_weight=weights)
                - roc_auc_score(y, anchor_scores, sample_weight=weights)
            )
            draws["AP"].append(
                average_precision_score(y, scores, sample_weight=weights)
                - average_precision_score(y, anchor_scores, sample_weight=weights)
            )
            draws["MCC_target_prior"].append(
                projected_mcc(y, calls, weights, cfg["target_prior"])
                - projected_mcc(y, anchor_calls, weights, cfg["target_prior"])
            )
        for metric, values in draws.items():
            low, median, high = np.quantile(values, [0.025, 0.5, 0.975])
            rows.append(
                {
                    "strategy": strategy,
                    "method": method,
                    "reference": cfg["anchor"],
                    "metric": metric,
                    "delta_p025": float(low),
                    "delta_p50": float(median),
                    "delta_p975": float(high),
                }
            )

result = pl.DataFrame(rows).sort("strategy", "method", "metric")
result.write_csv(out / "paired-group-bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "anchor": cfg["anchor"],
            "replicates": cfg["replicates"],
            "scope": (
                "Conditional paired group bootstrap of fixed OOF predictions and calls; "
                "excludes model-selection and refitting uncertainty."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(result.write_csv(), flush=True)
