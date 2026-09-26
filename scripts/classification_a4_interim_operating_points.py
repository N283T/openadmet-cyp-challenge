"""Infer approximate A4 scored-subset call rates from aggregate leaderboard metrics."""

import json
import os
from pathlib import Path

import polars as pl


def main():
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    board = pl.read_csv(inputs / config["leaderboard"])
    submission = pl.read_parquet(inputs / config["submission"])
    exact_calls = submission["CYP3A4_is_TDI"].sum()
    if submission.height != 750 or exact_calls != 249:
        raise ValueError("Unexpected submitted A4 call count")

    rows = []
    for user in config["users"]:
        match = board.filter(pl.col("Username") == user)
        if match.height != 1:
            raise ValueError(f"Missing or duplicate leaderboard row: {user}")
        item = match.row(0, named=True)
        precision = item["Precision"]
        recall = item["Recall"]
        accuracy = item["Accuracy"]
        denominator = 1 - 2 * recall + recall / precision
        prevalence = (1 - accuracy) / denominator
        call_fraction = prevalence * recall / precision
        rows.append(
            {
                "username": user,
                "mcc_rank": item["Rank"],
                "mcc": item["MCC"],
                "accuracy": accuracy,
                "precision": precision,
                "recall": recall,
                "f1": item["F1 Score"],
                "implied_scored_prevalence": prevalence,
                "implied_scored_call_fraction": call_fraction,
            }
        )
    pl.DataFrame(rows).write_csv(outputs / "operating_points.csv")
    (outputs / "summary.json").write_text(
        json.dumps(
            {
                "n283t_exact_test_positive_calls": exact_calls,
                "n283t_exact_test_call_fraction": exact_calls / submission.height,
                "caveat": (
                    "Leaderboard metrics are separate bootstrap means; derived prevalence "
                    "and scored-subset call fractions are approximations, not recovered "
                    "labels or exact submitted counts for other teams."
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
