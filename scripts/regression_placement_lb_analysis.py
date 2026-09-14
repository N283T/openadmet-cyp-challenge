"""Decompose the reflected all-endpoint placement leaderboard result."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read_table(filename: str) -> pl.DataFrame:
    path = inp / filename
    return pl.read_csv(path) if path.suffix == ".csv" else pl.read_parquet(path)


def participant_row(frame: pl.DataFrame, submitted: str) -> dict:
    rows = frame.filter(
        (pl.col("Username") == cfg["username"]) & (pl.col("Submitted") == submitted)
    )
    assert rows.height == 1, (submitted, rows.height)
    return rows.row(0, named=True)


result_rows = []
for endpoint in cfg["endpoints"]:
    old_board = read_table(endpoint["old"])
    new_board = read_table(endpoint["new"])
    old = participant_row(old_board, cfg["old_submitted"])
    new = participant_row(new_board, cfg["new_submitted"])
    score_column = endpoint["score_column"]
    result_rows.append(
        {
            "endpoint": endpoint["name"],
            "old_ST_RAE": float(old[score_column]),
            "new_ST_RAE": float(new[score_column]),
            "improvement": float(old[score_column] - new[score_column]),
            "old_rank": int(old["Rank"]),
            "new_rank": int(new["Rank"]),
            "rank_change": int(old["Rank"] - new["Rank"]),
            "leader_ST_RAE": float(new_board[score_column][0]),
            "remaining_gap_to_leader": float(new[score_column] - new_board[score_column][0]),
        }
    )

endpoint_results = pl.DataFrame(result_rows)
non_macro = endpoint_results.filter(pl.col("endpoint") != "MA")
total_endpoint_improvement = non_macro["improvement"].sum()
endpoint_results = endpoint_results.with_columns(
    pl.when(pl.col("endpoint") == "MA")
    .then(None)
    .otherwise(pl.col("improvement") / total_endpoint_improvement)
    .alias("share_of_endpoint_improvement")
)
endpoint_results.write_csv(out / "endpoint-results.csv")

scores = {
    row["endpoint"]: {"old": row["old_ST_RAE"], "new": row["new_ST_RAE"]}
    for row in non_macro.to_dicts()
}
counterfactual_rows = []
for scenario in cfg["counterfactuals"]:
    values = [
        scores[endpoint]["new" if endpoint in scenario["new_endpoints"] else "old"]
        for endpoint in cfg["endpoint_order"]
    ]
    counterfactual_rows.append(
        {
            "scenario": scenario["name"],
            "new_endpoints": "+".join(scenario["new_endpoints"]) or "none",
            "macro_ST_RAE": float(np.mean(values)),
            "improvement_vs_previous": float(
                np.mean([scores[endpoint]["old"] for endpoint in cfg["endpoint_order"]])
                - np.mean(values)
            ),
        }
    )
counterfactuals = pl.DataFrame(counterfactual_rows)
counterfactuals.write_csv(out / "counterfactuals.csv")

old_ma = participant_row(read_table(cfg["old_ma"]), cfg["old_submitted"])
new_ma = participant_row(read_table(cfg["new_ma"]), cfg["new_submitted"])
metric_rows = []
for metric in cfg["macro_metrics"]:
    old_value = float(old_ma[metric])
    new_value = float(new_ma[metric])
    metric_rows.append(
        {
            "metric": metric,
            "old": old_value,
            "new": new_value,
            "change": new_value - old_value,
        }
    )
pl.DataFrame(metric_rows).write_csv(out / "macro-metrics.csv")

old_predictions = read_table(cfg["old_predictions"]).sort("Molecule_Name")
new_predictions = read_table(cfg["new_predictions"]).sort("Molecule_Name")
assert old_predictions.select("Molecule_Name", "SMILES").equals(
    new_predictions.select("Molecule_Name", "SMILES")
)
transform_rows = []
for endpoint in cfg["endpoint_order"]:
    old = old_predictions[endpoint].to_numpy().astype(float)
    new = new_predictions[endpoint].to_numpy().astype(float)
    slope, intercept = np.polyfit(old, new, 1)
    reconstructed = slope * old + intercept
    transform_rows.append(
        {
            "endpoint": endpoint,
            "slope": float(slope),
            "intercept": float(intercept),
            "max_affine_error": float(np.max(np.abs(reconstructed - new))),
            "raw_mean": float(old.mean()),
            "placed_mean": float(new.mean()),
            "raw_sample_sd": float(old.std(ddof=1)),
            "placed_sample_sd": float(new.std(ddof=1)),
            "mean_absolute_change": float(np.abs(new - old).mean()),
            "rank_identical": bool(np.array_equal(rankdata(old), rankdata(new))),
        }
    )
transforms = pl.DataFrame(transform_rows)
transforms.write_csv(out / "transforms.csv")

observed_macro = float(
    endpoint_results.filter(pl.col("endpoint") == "MA")["new_ST_RAE"][0]
)
reconstructed_macro = float(non_macro["new_ST_RAE"].mean())
checks = {
    "endpoint_macro_reconstruction_error": abs(observed_macro - reconstructed_macro),
    "all_endpoint_scores_improved": bool((non_macro["improvement"] > 0).all()),
    "all_prediction_ranks_identical": bool(transforms["rank_identical"].all()),
    "maximum_affine_reconstruction_error": float(transforms["max_affine_error"].max()),
    "old_submission": cfg["old_submitted"],
    "new_submission": cfg["new_submitted"],
}
(out / "summary.json").write_text(
    json.dumps(
        {
            "checks": checks,
            "headline": {
                "old_macro_ST_RAE": float(old_ma["MA-ST-RAE"]),
                "new_macro_ST_RAE": float(new_ma["MA-ST-RAE"]),
                "old_rank": int(old_ma["Rank"]),
                "new_rank": int(new_ma["Rank"]),
                "total_improvement": float(old_ma["MA-ST-RAE"] - new_ma["MA-ST-RAE"]),
            },
            "interpretation_limits": [
                (
                    "Endpoint counterfactuals are exact under the leaderboard's "
                    "endpoint-wise macro mean."
                ),
                (
                    "The two submissions isolate affine placement because molecular "
                    "ranks are identical."
                ),
                (
                    "Hidden labels remain unavailable, so this result does not identify "
                    "the exact optimal future placement constants."
                ),
            ],
        },
        indent=2,
    )
    + "\n"
)
print(endpoint_results, flush=True)
print(counterfactuals, flush=True)
print(json.dumps(checks), flush=True)
