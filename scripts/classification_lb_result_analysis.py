"""Decompose the reflected diverse-model classification leaderboard gain."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoints = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]


def csv(name):
    return pl.read_csv(inp / cfg[name])


def parquet(name):
    return pl.read_parquet(inp / cfg[name]).sort("Molecule_Name")


def own_row(board, metric="MCC"):
    selected = board.filter(pl.col("Username") == cfg["our_username"])
    assert selected.height == 1
    row = selected.row(0, named=True)
    return {
        "rank": int(row["Rank"]),
        "score": float(row[metric]),
        "accuracy": float(row["Accuracy"] if metric == "MCC" else row["MA-Accuracy"]),
        "precision": float(row["Precision"] if metric == "MCC" else row["MA-Precision"]),
        "recall": float(row["Recall"] if metric == "MCC" else row["MA-Recall"]),
        "f1": float(row["F1 Score"] if metric == "MCC" else row["MA-F1 Score"]),
        "submitted": row["Submitted"],
    }


boards = {
    generation: {endpoint: csv(f"{generation}_{endpoint}") for endpoint in endpoints}
    for generation in ["old", "pre", "post"]
}
rows = []
priors = {}
for endpoint in endpoints:
    priors[endpoint] = cfg["leaderboard_priors"][endpoint]
    for generation in ["old", "pre", "post"]:
        row = own_row(boards[generation][endpoint])
        prevalence = priors[endpoint]
        implied_rate = row["recall"] * prevalence / row["precision"]
        true_positive_share = row["recall"] * prevalence
        rows.append(
            {
                "endpoint": endpoint,
                "generation": generation,
                **row,
                "hidden_prevalence": prevalence,
                "implied_positive_call_rate": implied_rate,
                "true_positive_share": true_positive_share,
                "false_positive_share": implied_rate - true_positive_share,
            }
        )
metrics = pl.DataFrame(rows)
metrics.write_csv(out / "endpoint-operating-points.csv")

deltas = []
for endpoint in endpoints:
    values = {
        row["generation"]: row
        for row in metrics.filter(pl.col("endpoint") == endpoint).iter_rows(named=True)
    }
    for comparison, left, right in [
        ("post_minus_pre", "post", "pre"),
        ("post_minus_old", "post", "old"),
    ]:
        deltas.append(
            {
                "endpoint": endpoint,
                "comparison": comparison,
                **{
                    f"delta_{field}": values[left][field] - values[right][field]
                    for field in [
                        "score",
                        "accuracy",
                        "precision",
                        "recall",
                        "f1",
                        "implied_positive_call_rate",
                        "true_positive_share",
                        "false_positive_share",
                    ]
                },
            }
        )
pl.DataFrame(deltas).write_csv(out / "endpoint-deltas.csv")

submissions = {
    generation: parquet(f"{generation}_submission") for generation in ["old", "pre", "post"]
}
identity = submissions["post"].select("SMILES", "Molecule_Name")
for submission in submissions.values():
    assert submission.select("SMILES", "Molecule_Name").equals(identity)
transition_rows = []
for endpoint in endpoints:
    for comparison, first, second in [
        ("pre_to_post", "pre", "post"),
        ("old_to_post", "old", "post"),
    ]:
        first_values = submissions[first][endpoint]
        second_values = submissions[second][endpoint]
        for first_call in [False, True]:
            for second_call in [False, True]:
                transition_rows.append(
                    {
                        "endpoint": endpoint,
                        "comparison": comparison,
                        "from_prediction": first_call,
                        "to_prediction": second_call,
                        "n": int(
                            ((first_values == first_call) & (second_values == second_call)).sum()
                        ),
                    }
                )
pl.DataFrame(transition_rows).write_csv(out / "prediction-transitions.csv")

scores = {
    endpoint: {
        generation: float(
            metrics.filter((pl.col("endpoint") == endpoint) & (pl.col("generation") == generation))[
                "score"
            ][0]
        )
        for generation in ["old", "pre", "post"]
    }
    for endpoint in endpoints
}
d6, a4 = endpoints
macro_rows = [
    {
        "combination": "pre_D6_pre_A4",
        "D6_MCC": scores[d6]["pre"],
        "A4_MCC": scores[a4]["pre"],
    },
    {
        "combination": "pre_D6_old_A4",
        "D6_MCC": scores[d6]["pre"],
        "A4_MCC": scores[a4]["old"],
    },
    {
        "combination": "post_D6_old_A4",
        "D6_MCC": scores[d6]["post"],
        "A4_MCC": scores[a4]["old"],
    },
    {
        "combination": "pre_D6_post_A4",
        "D6_MCC": scores[d6]["pre"],
        "A4_MCC": scores[a4]["post"],
    },
    {
        "combination": "post_D6_post_A4",
        "D6_MCC": scores[d6]["post"],
        "A4_MCC": scores[a4]["post"],
    },
]
macro = pl.DataFrame(macro_rows).with_columns(
    ((pl.col("D6_MCC") + pl.col("A4_MCC")) / 2).alias("macro_MCC")
)
macro.write_csv(out / "exact-endpoint-combinations.csv")

post_ma = own_row(csv("post_MA"), "MA-MCC")
pre_ma = own_row(csv("pre_MA"), "MA-MCC")
total_gain = post_ma["score"] - pre_ma["score"]
d6_gain = (scores[d6]["post"] - scores[d6]["pre"]) / 2
a4_gain = (scores[a4]["post"] - scores[a4]["pre"]) / 2
summary = {
    "pre_macro": pre_ma,
    "post_macro": post_ma,
    "macro_gain": total_gain,
    "D6_macro_contribution": d6_gain,
    "A4_macro_contribution": a4_gain,
    "D6_fraction_of_gain": d6_gain / total_gain,
    "A4_fraction_of_gain": a4_gain / total_gain,
    "tabicl_baseline": cfg["tabicl_baseline"],
    "margin_over_tabicl": post_ma["score"] - cfg["tabicl_baseline"],
    "leaderboard_priors": priors,
    "limitations": (
        "Leaderboard metrics are bootstrap means on the confidently scored subset; "
        "endpoint combinations are exact means of scored endpoint MCCs, while "
        "confusion shares and implied call rates are bootstrap-weighted estimates."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
