"""Review scored placements and package an OOF-guarded v9 correction."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = cfg["endpoint"]
keys = ["compound_id", "endpoint", "fold"]


def interval_score(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    y = frame["y_true"].to_numpy()
    low = frame["conf_low"].to_numpy()
    high = frame["conf_high"].to_numpy()
    error = np.maximum(low - prediction, 0) + np.maximum(prediction - high, 0)
    center = np.full_like(y, y.mean())
    denominator = (np.maximum(low - center, 0) + np.maximum(center - high, 0)).sum()
    if denominator <= 0:
        raise ValueError("Non-positive ST-RAE denominator")
    return float(error.sum() / denominator)


def read_method(filename: str, method: str, name: str) -> pl.DataFrame:
    return (
        pl.read_parquet(inp / filename)
        .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
        .select(*keys, "y_true", "conf_low", "conf_high", pl.col("prediction").alias(name))
        .sort(keys)
    )


old_oof = read_method("old-oof.parquet", cfg["old_oof_method"], "old")
new_oof = read_method("new-oof.parquet", cfg["new_oof_method"], "new")
oof = old_oof.join(
    new_oof.select(*keys, pl.col("new")), on=keys, validate="1:1", maintain_order="left"
)
if oof.height == 0 or oof["new"].null_count():
    raise ValueError("Incomplete OOF alignment")

alpha_rows = []
for alpha in cfg["delta_alpha_grid"]:
    prediction = oof["old"].to_numpy() + alpha * (oof["new"].to_numpy() - oof["old"].to_numpy())
    alpha_rows.append(
        {
            "alpha": alpha,
            "fold": -1,
            "ST_RAE": interval_score(oof, prediction),
            "MAE": float(np.abs(oof["y_true"].to_numpy() - prediction).mean()),
            "Spearman": float(spearmanr(oof["y_true"].to_numpy(), prediction).statistic),
        }
    )
    for fold in sorted(oof["fold"].unique().to_list()):
        mask = oof["fold"].to_numpy() == fold
        alpha_rows.append(
            {
                "alpha": alpha,
                "fold": fold,
                "ST_RAE": interval_score(oof.filter(pl.col("fold") == fold), prediction[mask]),
                "MAE": float(np.abs(oof["y_true"].to_numpy()[mask] - prediction[mask]).mean()),
                "Spearman": float(
                    spearmanr(oof["y_true"].to_numpy()[mask], prediction[mask]).statistic
                ),
            }
        )
alpha_metrics = pl.DataFrame(alpha_rows).sort("fold", "alpha")
alpha_metrics.write_csv(out / "d6-delta-alpha-metrics.csv")

baseline_folds = alpha_metrics.filter((pl.col("alpha") == 1.0) & (pl.col("fold") >= 0)).select(
    "fold", pl.col("ST_RAE").alias("baseline")
)
eligible = (
    alpha_metrics.filter((pl.col("fold") >= 0) & (pl.col("alpha") <= cfg["selection_alpha_max"]))
    .join(baseline_folds, on="fold", validate="m:1")
    .group_by("alpha")
    .agg(
        (pl.col("ST_RAE") - pl.col("baseline")).max().alias("worst_fold_delta"),
    )
    .filter(pl.col("worst_fold_delta") <= cfg["fold_guardrail"])
    .join(alpha_metrics.filter(pl.col("fold") == -1).select("alpha", "ST_RAE"), on="alpha")
    .sort("ST_RAE", "alpha")
)
if eligible.height == 0:
    selected_alpha = 1.0
    selection_reason = "fallback_alpha_1_no_guardrail_eligible_variant"
else:
    selected_alpha = float(eligible["alpha"][0])
    selection_reason = "minimum_pooled_oof_with_fold_guardrail_and_alpha_cap"


def read_predictions(filename: str) -> pl.DataFrame:
    path = inp / filename
    return pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(path)


old_raw = read_predictions("old-raw.csv").sort("Molecule_Name")
new_raw = read_predictions("new-raw.csv").sort("Molecule_Name")
current = read_predictions("current-placed.csv").sort("Molecule_Name")
reference = pl.read_csv(inp / "reference.csv")
for frame in [new_raw, current]:
    if not frame["Molecule_Name"].equals(old_raw["Molecule_Name"]):
        raise ValueError("Prediction identity mismatch")

corrected = new_raw
raw_d6 = old_raw[endpoint].to_numpy() + selected_alpha * (
    new_raw[endpoint].to_numpy() - old_raw[endpoint].to_numpy()
)
corrected = corrected.with_columns(pl.Series(endpoint, raw_d6))

placement_rows = []
for column, target in cfg["placements"].items():
    values = corrected[column].to_numpy()
    adjusted = (values - values.mean()) / values.std(ddof=1) * target["sample_sd"] + target["mean"]
    corrected = corrected.with_columns(pl.Series(column, adjusted))
    placement_rows.append(
        {
            "endpoint": column,
            "raw_mean": float(values.mean()),
            "raw_sample_sd": float(values.std(ddof=1)),
            "target_mean": target["mean"],
            "target_sample_sd": target["sample_sd"],
        }
    )

corrected.write_csv(out / "candidate-placement-v9.csv")
corrected.write_parquet(out / "candidate-placement-v9.parquet")
report = validate(out / "candidate-placement-v9.csv", "regression", inp / "reference.csv")
if not report["valid"]:
    raise ValueError(report)

d6_oof_pearson = float(np.corrcoef(oof["y_true"].to_numpy(), oof["new"].to_numpy())[0, 1])
d6_correlation_sd = d6_oof_pearson * cfg["d6_inferred_blind_sd"]
d6_sd_choices = cfg["d6_placement_sensitivity"]["sample_sds"]
d6_adaptive_sd = min(d6_sd_choices, key=lambda value: abs(value - d6_correlation_sd))
adaptive = corrected.with_columns(
    (
        cfg["placements"][endpoint]["mean"]
        + (pl.col(endpoint) - pl.col(endpoint).mean())
        / pl.col(endpoint).std(ddof=1)
        * d6_adaptive_sd
    ).alias(endpoint)
)
adaptive.write_csv(out / "candidate-placement-v9-adaptive-d6.csv")
adaptive.write_parquet(out / "candidate-placement-v9-adaptive-d6.parquet")
adaptive_report = validate(
    out / "candidate-placement-v9-adaptive-d6.csv", "regression", inp / "reference.csv"
)
if not adaptive_report["valid"]:
    raise ValueError(adaptive_report)

comparison_rows = []
for column in cfg["placements"]:
    before = current[column].to_numpy()
    after = corrected[column].to_numpy()
    comparison_rows.append(
        {
            "endpoint": column,
            "mean_absolute_change": float(np.abs(after - before).mean()),
            "maximum_absolute_change": float(np.abs(after - before).max()),
            "spearman": float(spearmanr(before, after).statistic),
            "old_mean": float(before.mean()),
            "new_mean": float(after.mean()),
            "old_sample_sd": float(before.std(ddof=1)),
            "new_sample_sd": float(after.std(ddof=1)),
        }
    )
comparison = pl.DataFrame(comparison_rows)
comparison.write_csv(out / "comparison-to-v8.csv")
pl.DataFrame(placement_rows).write_csv(out / "placement-moments.csv")

sensitivity_rows = []
z = (raw_d6 - raw_d6.mean()) / raw_d6.std(ddof=1)
for mean in cfg["d6_placement_sensitivity"]["means"]:
    for sample_sd in cfg["d6_placement_sensitivity"]["sample_sds"]:
        values = mean + sample_sd * z
        sensitivity_rows.append(
            {
                "mean": mean,
                "sample_sd": sample_sd,
                "mean_absolute_change_vs_current": float(
                    np.abs(values - current[endpoint].to_numpy()).mean()
                ),
                "minimum": float(values.min()),
                "maximum": float(values.max()),
            }
        )
pl.DataFrame(sensitivity_rows).write_csv(out / "d6-placement-sensitivity.csv")

history_rows = []
prior = None
for item in cfg["leaderboard_history"]:
    frame = read_predictions(item["file"]).sort("Molecule_Name")
    for column in cfg["placements"]:
        values = frame[column].to_numpy()
        key = column.split("_")[0].replace("CYP", "")
        history_rows.append(
            {
                "submission": item["name"],
                "endpoint": column,
                "ST_RAE": item["scores"][key],
                "mean": float(values.mean()),
                "sample_sd": float(values.std(ddof=1)),
                "spearman_to_previous": (
                    None if prior is None else float(spearmanr(prior[column], values).statistic)
                ),
            }
        )
    prior = {column: frame[column].to_numpy() for column in cfg["placements"]}
pl.DataFrame(history_rows).write_csv(out / "leaderboard-history.csv")

summary = {
    "selected_alpha": selected_alpha,
    "selection_reason": selection_reason,
    "eligible_alphas": eligible.to_dicts(),
    "placement": placement_rows,
    "comparison_to_v8": comparison_rows,
    "validation": report,
    "adaptive_d6": {
        "oof_pearson": d6_oof_pearson,
        "inferred_blind_sd": cfg["d6_inferred_blind_sd"],
        "correlation_scaled_sd": d6_correlation_sd,
        "selected_sample_sd": d6_adaptive_sd,
        "candidate": "candidate-placement-v9-adaptive-d6.csv",
        "validation": adaptive_report,
        "status": "adventurous_submission_candidate",
    },
    "leaderboard_scores_are_observed_not_estimated": True,
    "test_labels_used": False,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
