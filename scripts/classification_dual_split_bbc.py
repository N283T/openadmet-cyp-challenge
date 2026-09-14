"""Group-bootstrap candidate-selection audit for the compact dual-split screen."""

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["oof"])


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(positions))))
        order = np.argsort(-scores[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def expected_mcc(
    y: np.ndarray, scores: np.ndarray, folds: np.ndarray, rate: float, prevalence: float
) -> float:
    calls = fixed_rate_calls(scores, folds, rate)
    positive = y == 1
    if not positive.any() or positive.all():
        return float("nan")
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


rng = np.random.default_rng(cfg["seed"])
replicate_rows = []
summary_rows = []
frequency_rows = []
for endpoint in cfg["endpoints"]:
    settings = cfg["endpoint_settings"][endpoint]
    for strategy in cfg["strategies"]:
        table = oof.filter((pl.col("endpoint") == endpoint) & (pl.col("strategy") == strategy))
        methods = sorted(table["method"].unique().to_list())
        base = table.filter(pl.col("method") == methods[0]).sort("compound_id")
        y = base["y_true"].to_numpy().astype(int)
        folds = base["fold"].to_numpy().astype(int)
        groups = base["group"].cast(pl.String).to_numpy()
        scores = {}
        for method in methods:
            selected = table.filter(pl.col("method") == method).sort("compound_id")
            if selected["compound_id"].to_list() != base["compound_id"].to_list():
                raise ValueError(f"Prediction keys differ for {method}")
            scores[method] = selected["prediction"].to_numpy()

        full_scores = {
            method: expected_mcc(
                y,
                prediction,
                folds,
                settings["call_rate"],
                settings["target_prior"],
            )
            for method, prediction in scores.items()
        }
        full_winner = max(full_scores, key=full_scores.get)
        group_values = np.unique(groups)
        group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
        selection_counts = Counter()
        oob_scores = []
        apparent_scores = []
        for replicate in range(cfg["replicates"]):
            sampled_groups = rng.choice(group_values, size=len(group_values), replace=True)
            sampled_set = set(sampled_groups)
            selected_rows = np.concatenate([group_rows[group] for group in sampled_groups])
            oob_groups = [group for group in group_values if group not in sampled_set]
            if not oob_groups:
                continue
            oob_rows = np.concatenate([group_rows[group] for group in oob_groups])
            candidate_scores = {
                method: expected_mcc(
                    y[selected_rows],
                    prediction[selected_rows],
                    folds[selected_rows],
                    settings["call_rate"],
                    settings["target_prior"],
                )
                for method, prediction in scores.items()
            }
            winner = max(candidate_scores, key=candidate_scores.get)
            oob_score = expected_mcc(
                y[oob_rows],
                scores[winner][oob_rows],
                folds[oob_rows],
                settings["call_rate"],
                settings["target_prior"],
            )
            if not np.isfinite(oob_score):
                continue
            selection_counts[winner] += 1
            apparent_scores.append(candidate_scores[winner])
            oob_scores.append(oob_score)
            replicate_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "replicate": replicate,
                    "winner": winner,
                    "apparent_selection_score": candidate_scores[winner],
                    "oob_score": oob_score,
                    "oob_groups": len(oob_groups),
                    "oob_rows": len(oob_rows),
                }
            )
        oob_values = np.asarray(oob_scores)
        apparent_values = np.asarray(apparent_scores)
        summary_rows.append(
            {
                "endpoint": endpoint,
                "strategy": strategy,
                "full_oof_winner": full_winner,
                "full_oof_winner_score": full_scores[full_winner],
                "valid_replicates": len(oob_values),
                "bbc_oob_mean": float(oob_values.mean()),
                "bbc_oob_p025": float(np.quantile(oob_values, 0.025)),
                "bbc_oob_p975": float(np.quantile(oob_values, 0.975)),
                "mean_apparent_optimism": float((apparent_values - oob_values).mean()),
                "winner_selection_frequency": selection_counts[full_winner] / len(oob_values),
            }
        )
        for method in methods:
            frequency_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "method": method,
                    "full_oof_score": full_scores[method],
                    "bootstrap_wins": selection_counts[method],
                    "bootstrap_win_fraction": selection_counts[method] / len(oob_values),
                }
            )

summary = pl.DataFrame(summary_rows).sort("endpoint", "strategy")
summary.write_csv(out / "summary.csv")
pl.DataFrame(frequency_rows).sort(
    "endpoint", "strategy", "bootstrap_win_fraction", descending=[False, False, True]
).write_csv(out / "selection-frequency.csv")
pl.DataFrame(replicate_rows).write_parquet(out / "bootstrap-replicates.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "method": "group bootstrap with out-of-bag evaluation after candidate selection",
            "replicates": cfg["replicates"],
            "seed": cfg["seed"],
            "scope": "Selection-bias diagnostic conditional on saved OOF predictions",
        },
        indent=2,
    )
    + "\n"
)
print(summary.write_csv(), flush=True)
