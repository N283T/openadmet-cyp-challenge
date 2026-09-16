"""Analyze the reflected D6 mean-block plus A4 equal-seven submission."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoints = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]


def csv(name: str) -> pl.DataFrame:
    return pl.read_csv(inp / cfg[name])


def own_row(board: pl.DataFrame, endpoint: str) -> dict:
    row = board.filter(pl.col("Username") == cfg["username"]).row(0, named=True)
    prefix = "MA-" if endpoint == "MA" else ""
    return {
        "rank": int(row["Rank"]),
        "MCC": float(row[f"{prefix}MCC"]),
        "MCC_std": float(row[f"{prefix}MCC std"]),
        "accuracy": float(row[f"{prefix}Accuracy"]),
        "precision": float(row[f"{prefix}Precision"]),
        "recall": float(row[f"{prefix}Recall"]),
        "f1": float(row[f"{prefix}F1 Score"]),
        "submitted": row["Submitted"],
    }


boards = {
    generation: {
        endpoint: csv(f"{generation}_{endpoint}") for endpoint in ["MA", *endpoints]
    }
    for generation in ["pre", "post"]
}
metric_rows = []
delta_rows = []
for endpoint in ["MA", *endpoints]:
    values = {}
    for generation in ["pre", "post"]:
        values[generation] = own_row(boards[generation][endpoint], endpoint)
        metric_rows.append(
            {"endpoint": endpoint, "generation": generation, **values[generation]}
        )
    delta_rows.append(
        {
            "endpoint": endpoint,
            "rank_change": values["post"]["rank"] - values["pre"]["rank"],
            **{
                f"delta_{field}": values["post"][field] - values["pre"][field]
                for field in ["MCC", "accuracy", "precision", "recall", "f1"]
            },
        }
    )
pl.DataFrame(metric_rows).write_csv(out / "leaderboard-metrics.csv")
pl.DataFrame(delta_rows).write_csv(out / "leaderboard-deltas.csv")

pre_submission = pl.read_parquet(inp / cfg["pre_submission"]).sort("Molecule_Name")
post_submission = pl.read_parquet(inp / cfg["post_submission"]).sort("Molecule_Name")
if not pre_submission.select("SMILES", "Molecule_Name").equals(
    post_submission.select("SMILES", "Molecule_Name")
):
    raise ValueError("Submission identities differ")
transition_rows = []
for endpoint in endpoints:
    before, after = pre_submission[endpoint], post_submission[endpoint]
    transition_rows.append(
        {
            "endpoint": endpoint,
            "pre_positive": int(before.sum()),
            "post_positive": int(after.sum()),
            "both_positive": int((before & after).sum()),
            "pre_only": int((before & ~after).sum()),
            "post_only": int((~before & after).sum()),
        }
    )
pl.DataFrame(transition_rows).write_csv(out / "prediction-transitions.csv")

post_scores = {
    endpoint: own_row(boards["post"][endpoint], endpoint)["MCC"]
    for endpoint in endpoints
}
pre_scores = {
    endpoint: own_row(boards["pre"][endpoint], endpoint)["MCC"]
    for endpoint in endpoints
}
macro_gain = (sum(post_scores.values()) - sum(pre_scores.values())) / 2
contributions = {
    endpoint: (post_scores[endpoint] - pre_scores[endpoint]) / 2 for endpoint in endpoints
}

macro_board = boards["post"]["MA"]
a4_board = boards["post"]["CYP3A4_is_TDI"]
macro_us = own_row(macro_board, "MA")
a4_us = own_row(a4_board, "CYP3A4_is_TDI")
rank1_macro = float(macro_board.filter(pl.col("Rank") == 1)["MA-MCC"][0])
rank2_macro = float(macro_board.filter(pl.col("Rank") == 2)["MA-MCC"][0])
rank4_macro = float(macro_board.filter(pl.col("Rank") == 4)["MA-MCC"][0])
rank2_a4 = float(a4_board.filter(pl.col("Rank") == 2)["MCC"][0])

a4_bootstrap = csv("a4_bootstrap").filter(
    pl.col("model") == "exhaustive_mcc::components_7"
).row(0, named=True)
d6_metrics = csv("d6_metrics")
d6_pca = d6_metrics.filter(pl.col("model") == "crossfit::pca256").row(0, named=True)
d6_mean = d6_metrics.filter(pl.col("model") == "crossfit::mean_block").row(0, named=True)

summary = {
    "reflected_submission": own_row(boards["post"]["MA"], "MA")["submitted"],
    "macro": {
        "pre": own_row(boards["pre"]["MA"], "MA"),
        "post": macro_us,
        "gain": macro_gain,
        "margin_to_rank_2": macro_us["MCC"] - rank2_macro,
        "margin_over_rank_4": macro_us["MCC"] - rank4_macro,
    },
    "endpoint_MCC": {
        endpoint: {"pre": pre_scores[endpoint], "post": post_scores[endpoint]}
        for endpoint in endpoints
    },
    "macro_gain_contribution": contributions,
    "macro_gain_fraction": {
        endpoint: value / macro_gain for endpoint, value in contributions.items()
    },
    "A4_rank_1_margin": a4_us["MCC"] - rank2_a4,
    "D6_MCC_needed_for_current_macro_rank_1_score_if_A4_fixed": (
        2 * rank1_macro - post_scores["CYP3A4_is_TDI"]
    ),
    "D6_additional_MCC_needed": 2 * rank1_macro - sum(post_scores.values()),
    "CV_context": {
        "A4_equal7_target_MCC_delta": float(a4_bootstrap["mean_MCC_delta"]),
        "A4_equal7_group_CI": [float(a4_bootstrap["ci025"]), float(a4_bootstrap["ci975"])],
        "D6_PCA_to_mean_block_target_MCC_delta": float(
            d6_mean["MCC_target_prior"] - d6_pca["MCC_target_prior"]
        ),
        "D6_PCA_to_mean_block_AUROC_delta": float(d6_mean["AUROC"] - d6_pca["AUROC"]),
    },
    "interpretation_limits": [
        (
            "Both endpoints retained identical positive counts, so the gains reflect "
            "changed rankings/calls rather than call-rate changes."
        ),
        (
            "The A4 submission changed both ensemble composition and weighting; one "
            "leaderboard result cannot attribute the gain specifically to raw embeddings."
        ),
        (
            "The reported leaderboard standard deviations are marginal bootstrap "
            "dispersions, not paired uncertainty for score differences or ranks."
        ),
        (
            "The A4 leaderboard gain is much larger than its grouped-CV estimate, so "
            "winner's curse or evaluation-set alignment remains plausible."
        ),
    ],
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
