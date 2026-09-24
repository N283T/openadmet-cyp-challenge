"""Evaluate all-molecule versus labelled-only D6 decision boundaries."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.classification_evaluation import (
    evaluate_rates,
    fixed_rate_calls,
    prepare_predictions,
    projected_mcc,
)

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
prediction_files = cfg["predictions"]
if isinstance(prediction_files, str):
    prediction_files = [prediction_files]
bank = pl.concat([pl.read_parquet(inputs / filename) for filename in prediction_files])
assignments = pl.read_parquet(inputs / cfg["assignments"])
if not set(cfg["models"]).issubset(bank["model"].unique().to_list()):
    raise ValueError("One or more configured models are missing")
bank = bank.filter(pl.col("model").is_in(cfg["models"]))
rows = []
fixed_rows = []
for strategy in cfg["strategies"]:
    full = bank.filter(pl.col("split_id") == strategy)
    for scope, predictions in (
        ("full_pool", full),
        ("labelled_only_diagnostic", full.filter(pl.col("y_true").is_not_null())),
    ):
        data = prepare_predictions(predictions, assignments, cfg["models"], strategy)
        reference = data.filter(pl.col("model") == cfg["models"][0]).sort("compound_id")
        folds = reference["fold"].to_numpy()
        groups = reference["group"].to_numpy()
        y_all = reference["y_true"].to_numpy()
        scored = ~np.isnan(y_all.astype(float))
        truth = y_all[scored].astype(int)
        fixed = {
            model: fixed_rate_calls(
                data.filter(pl.col("model") == model).sort("compound_id")["score"].to_numpy(),
                folds,
                cfg["historical_rate"],
            )
            for model in cfg["models"]
        }
        unique_groups = np.unique(groups[scored])
        group_rows = {group: np.flatnonzero(scored & (groups == group)) for group in unique_groups}
        for prior in cfg["target_priors"]:
            fixed_scores = {
                model: projected_mcc(truth, calls[scored], prior) for model, calls in fixed.items()
            }
            rng = np.random.default_rng(cfg["bootstrap_seed"])
            sampled_deltas = []
            for _ in range(cfg["bootstrap_replicates"]):
                sampled_groups = rng.choice(unique_groups, len(unique_groups), replace=True)
                sampled_rows = np.concatenate([group_rows[group] for group in sampled_groups])
                sampled_truth = y_all[sampled_rows].astype(int)
                delta = projected_mcc(
                    sampled_truth, fixed[cfg["models"][1]][sampled_rows], prior
                ) - projected_mcc(sampled_truth, fixed[cfg["models"][0]][sampled_rows], prior)
                if np.isfinite(delta):
                    sampled_deltas.append(delta)
            if not sampled_deltas:
                raise ValueError("No valid fixed-rate bootstrap samples")
            fixed_rows.append(
                {
                    "strategy": strategy,
                    "scope": scope,
                    "prior": prior,
                    "rate": cfg["historical_rate"],
                    "pool_rows": len(folds),
                    "scored_rows": int(scored.sum()),
                    "control_mcc": fixed_scores[cfg["models"][0]],
                    "candidate_mcc": fixed_scores[cfg["models"][1]],
                    "candidate_minus_control": fixed_scores[cfg["models"][1]]
                    - fixed_scores[cfg["models"][0]],
                    "delta_ci025": float(np.quantile(sampled_deltas, 0.025)),
                    "delta_ci975": float(np.quantile(sampled_deltas, 0.975)),
                    "control_calls_all": int(fixed[cfg["models"][0]].sum()),
                    "candidate_calls_all": int(fixed[cfg["models"][1]].sum()),
                    "control_calls_scored": int(fixed[cfg["models"][0]][scored].sum()),
                    "candidate_calls_scored": int(fixed[cfg["models"][1]][scored].sum()),
                }
            )
            summary, curve, decisions, calls = evaluate_rates(
                data,
                cfg["models"],
                cfg["rates"],
                prior,
                cfg["expected_pool_rows"],
                cfg["bootstrap_replicates"],
                cfg["bootstrap_seed"],
                allow_partial_pool=scope == "labelled_only_diagnostic",
            )
            tag = f"{strategy}-{scope}-prior{prior:.5f}"
            curve.write_csv(outputs / f"{tag}-curve.csv")
            decisions.write_csv(outputs / f"{tag}-choices.csv")
            calls.write_parquet(outputs / f"{tag}-calls.parquet")
            (outputs / f"{tag}-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            delta = summary["paired_delta_candidate_minus_control"]
            rows.append(
                {
                    "strategy": strategy,
                    "scope": scope,
                    "prior": prior,
                    "pool_rows": summary["pool_rows"],
                    "scored_rows": summary["scored_rows"],
                    "candidate_minus_control": delta["point"],
                    "delta_ci025": delta["ci025"],
                    "delta_ci975": delta["ci975"],
                    "rate_choices_control": json.dumps(
                        summary["nested_evaluation"][0]["selected_rates"]
                    ),
                    "rate_choices_candidate": json.dumps(
                        summary["nested_evaluation"][1]["selected_rates"]
                    ),
                    "model_choices": json.dumps(summary["selector_choices"]),
                }
            )
pl.DataFrame(rows).write_csv(outputs / "comparison.csv")
pl.DataFrame(fixed_rows).write_csv(outputs / "fixed-rate.csv")
print(
    pl.DataFrame(rows).select(
        "strategy",
        "scope",
        "prior",
        "pool_rows",
        "scored_rows",
        "candidate_minus_control",
        "delta_ci025",
        "delta_ci975",
    ),
    flush=True,
)
