"""Paired-bootstrap analysis of the structural-alert ablation."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["oof"])

pairs = (
    ("morgan_balanced", "morgan_alerts_balanced"),
    ("fp_mordred_balanced", "fp_mordred_alerts_balanced"),
)


def score(y, probability, threshold):
    return {
        "MCC": float(matthews_corrcoef(y, probability >= threshold)),
        "AUROC": float(roc_auc_score(y, probability)),
        "AP": float(average_precision_score(y, probability)),
    }


rows = []
bootstrap_deltas = {}
for pair_index, (baseline, augmented) in enumerate(pairs):
    for policy in ("fixed_05", "inner_mcc"):
        for endpoint_index, endpoint in enumerate(cfg["endpoints"]):
            base = (
                oof.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == baseline))
                .sort("compound_id")
            )
            alert = (
                oof.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == augmented))
                .sort("compound_id")
            )
            assert base["compound_id"].equals(alert["compound_id"])
            assert base["y_true"].equals(alert["y_true"])
            y = base["y_true"].to_numpy().astype(int)
            base_probability = base["prediction"].to_numpy()
            alert_probability = alert["prediction"].to_numpy()
            base_threshold = 0.5 if policy == "fixed_05" else base["threshold"].to_numpy()
            alert_threshold = 0.5 if policy == "fixed_05" else alert["threshold"].to_numpy()
            base_score = score(y, base_probability, base_threshold)
            alert_score = score(y, alert_probability, alert_threshold)

            rng = np.random.default_rng(20260913 + pair_index * 100 + endpoint_index * 10)
            deltas = np.empty(cfg["bootstrap_replicates"], dtype=float)
            for replicate in range(cfg["bootstrap_replicates"]):
                sample = rng.integers(0, len(y), len(y))
                base_cut = base_threshold if np.isscalar(base_threshold) else base_threshold[sample]
                alert_cut = (
                    alert_threshold if np.isscalar(alert_threshold) else alert_threshold[sample]
                )
                deltas[replicate] = matthews_corrcoef(
                    y[sample], alert_probability[sample] >= alert_cut
                ) - matthews_corrcoef(y[sample], base_probability[sample] >= base_cut)
            bootstrap_deltas[(baseline, policy, endpoint)] = deltas
            rows.append(
                {
                    "endpoint": endpoint,
                    "policy": policy,
                    "baseline": baseline,
                    "augmented": augmented,
                    "baseline_MCC": base_score["MCC"],
                    "augmented_MCC": alert_score["MCC"],
                    "delta_MCC": alert_score["MCC"] - base_score["MCC"],
                    "delta_MCC_ci_low": float(np.quantile(deltas, 0.025)),
                    "delta_MCC_ci_high": float(np.quantile(deltas, 0.975)),
                    "probability_delta_MCC_gt_zero": float((deltas > 0).mean()),
                    "baseline_AUROC": base_score["AUROC"],
                    "augmented_AUROC": alert_score["AUROC"],
                    "delta_AUROC": alert_score["AUROC"] - base_score["AUROC"],
                    "baseline_AP": base_score["AP"],
                    "augmented_AP": alert_score["AP"],
                    "delta_AP": alert_score["AP"] - base_score["AP"],
                }
            )

for baseline, augmented in pairs:
    for policy in ("fixed_05", "inner_mcc"):
        endpoint_rows = [
            row
            for row in rows
            if row["baseline"] == baseline and row["policy"] == policy
        ]
        macro_deltas = np.mean(
            [
                bootstrap_deltas[(baseline, policy, endpoint)]
                for endpoint in cfg["endpoints"]
            ],
            axis=0,
        )
        rows.append(
            {
                "endpoint": "macro",
                "policy": policy,
                "baseline": baseline,
                "augmented": augmented,
                "baseline_MCC": float(np.mean([row["baseline_MCC"] for row in endpoint_rows])),
                "augmented_MCC": float(
                    np.mean([row["augmented_MCC"] for row in endpoint_rows])
                ),
                "delta_MCC": float(np.mean([row["delta_MCC"] for row in endpoint_rows])),
                "delta_MCC_ci_low": float(np.quantile(macro_deltas, 0.025)),
                "delta_MCC_ci_high": float(np.quantile(macro_deltas, 0.975)),
                "probability_delta_MCC_gt_zero": float((macro_deltas > 0).mean()),
                "baseline_AUROC": float(
                    np.mean([row["baseline_AUROC"] for row in endpoint_rows])
                ),
                "augmented_AUROC": float(
                    np.mean([row["augmented_AUROC"] for row in endpoint_rows])
                ),
                "delta_AUROC": float(np.mean([row["delta_AUROC"] for row in endpoint_rows])),
                "baseline_AP": float(np.mean([row["baseline_AP"] for row in endpoint_rows])),
                "augmented_AP": float(np.mean([row["augmented_AP"] for row in endpoint_rows])),
                "delta_AP": float(np.mean([row["delta_AP"] for row in endpoint_rows])),
            }
        )

result = pl.DataFrame(rows).sort("baseline", "policy", "endpoint")
result.write_csv(out / "paired-comparison.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "source_run": cfg["source_run"],
            "bootstrap_replicates": cfg["bootstrap_replicates"],
            "resampling_unit": "compound within endpoint",
            "interval": "paired percentile 95% interval",
        },
        indent=2,
    )
    + "\n"
)
print(result.write_csv(), flush=True)
