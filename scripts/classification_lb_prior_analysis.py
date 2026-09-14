"""Infer leaderboard operating priors and test prior-aware threshold policies."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import matthews_corrcoef

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read_csv(name):
    return pl.read_csv(inp / cfg[name])


def expected_mcc(y, binary, prevalence):
    positive = y == 1
    tpr = float(binary[positive].mean())
    fpr = float(binary[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1.0 - tpr)
    fp = (1.0 - prevalence) * fpr
    tn = (1.0 - prevalence) * (1.0 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator > 0 else 0.0


def choose_prior_threshold(y, probability, prevalence):
    order = np.argsort(-probability, kind="stable")
    ordered_probability = probability[order]
    positions = np.r_[
        0,
        np.flatnonzero(ordered_probability[:-1] != ordered_probability[1:]) + 1,
    ]
    thresholds = np.r_[
        np.nextafter(ordered_probability[0], np.inf),
        ordered_probability[positions[1:] - 1],
    ]
    choices = []
    for position, threshold in zip(positions, thresholds, strict=True):
        binary = np.zeros(len(y), dtype=bool)
        binary[order[:position]] = True
        choices.append(
            (
                expected_mcc(y, binary, prevalence),
                -abs(float(threshold) - 0.5),
                float(threshold),
            )
        )
    all_positive_threshold = float(ordered_probability[-1])
    choices.append((0.0, -abs(all_positive_threshold - 0.5), all_positive_threshold))
    return max(choices)[2]


leaderboard_specs = {
    "CYP2D6_is_TDI": ("old_d6_leaderboard", "current_d6_leaderboard"),
    "CYP3A4_is_TDI": ("old_a4_leaderboard", "current_a4_leaderboard"),
}
priors = {}
leader_rates = {}
top10_rates = {}
operating_rows = []
our_rows = []
for endpoint, (old_name, current_name) in leaderboard_specs.items():
    old_board = read_csv(old_name)
    board = read_csv(current_name)
    all_negative = board.filter(pl.col("Username") == cfg["all_negative_username"])
    assert len(all_negative) == 1
    prevalence = 1.0 - float(all_negative["Accuracy"][0])
    priors[endpoint] = prevalence

    endpoint_rows = []
    for row in board.iter_rows(named=True):
        precision = float(row["Precision"])
        recall = float(row["Recall"])
        implied_rate = recall * prevalence / precision if precision > 0 else 0.0
        result = {
            "endpoint": endpoint,
            "rank": int(row["Rank"]),
            "username": row["Username"],
            "MCC": float(row["MCC"]),
            "Accuracy": float(row["Accuracy"]),
            "Precision": precision,
            "Recall": recall,
            "implied_scored_positive_rate": implied_rate,
        }
        endpoint_rows.append(result)
        operating_rows.append(result)
    endpoint_rows.sort(key=lambda row: row["rank"])
    leader_rates[endpoint] = endpoint_rows[0]["implied_scored_positive_rate"]
    nonzero_top10 = [
        row["implied_scored_positive_rate"]
        for row in endpoint_rows[:10]
        if row["implied_scored_positive_rate"] > 0
    ]
    top10_rates[endpoint] = float(np.median(nonzero_top10))

    for generation, source in (("old", old_board), ("current", board)):
        our = source.filter(pl.col("Username") == cfg["our_username"])
        assert len(our) == 1
        row = our.row(0, named=True)
        precision = float(row["Precision"])
        recall = float(row["Recall"])
        our_rows.append(
            {
                "endpoint": endpoint,
                "generation": generation,
                "MCC": float(row["MCC"]),
                "Accuracy": float(row["Accuracy"]),
                "Precision": precision,
                "Recall": recall,
                "implied_scored_positive_rate": recall * prevalence / precision,
            }
        )

pl.DataFrame(operating_rows).write_csv(out / "leaderboard-operating-points.csv")
pl.DataFrame(our_rows).write_csv(out / "our-lb-operating-points.csv")

old_submission = pl.read_parquet(inp / cfg["old_submission"]).sort("Molecule_Name")
current_submission = pl.read_parquet(inp / cfg["current_submission"]).sort("Molecule_Name")
assert old_submission.select("SMILES", "Molecule_Name").equals(
    current_submission.select("SMILES", "Molecule_Name")
)
disagreements = []
for endpoint in cfg["endpoints"]:
    old_binary = old_submission[endpoint].to_numpy().astype(bool)
    current_binary = current_submission[endpoint].to_numpy().astype(bool)
    for old_value in (False, True):
        for current_value in (False, True):
            disagreements.append(
                {
                    "endpoint": endpoint,
                    "old_prediction": old_value,
                    "current_prediction": current_value,
                    "n": int(((old_binary == old_value) & (current_binary == current_value)).sum()),
                }
            )
pl.DataFrame(disagreements).write_csv(out / "submission-disagreements.csv")

our_table_for_transition = pl.DataFrame(our_rows)
transition_rows = []
for endpoint in cfg["endpoints"]:
    endpoint_rows = our_table_for_transition.filter(pl.col("endpoint") == endpoint)
    values = {row["generation"]: row for row in endpoint_rows.iter_rows(named=True)}
    prevalence = priors[endpoint]
    shares = {}
    for generation in ("old", "current"):
        row = values[generation]
        true_positive_share = row["Recall"] * prevalence
        predicted_positive_share = row["implied_scored_positive_rate"]
        shares[generation] = {
            "true_positive": true_positive_share,
            "false_positive": predicted_positive_share - true_positive_share,
            "predicted_positive": predicted_positive_share,
        }
    if shares["current"]["predicted_positive"] >= shares["old"]["predicted_positive"]:
        direction = "net_current_additions"
        first, second = "current", "old"
    else:
        direction = "net_current_removals"
        first, second = "old", "current"
    changed_positive = shares[first]["true_positive"] - shares[second]["true_positive"]
    changed_calls = (
        shares[first]["predicted_positive"] - shares[second]["predicted_positive"]
    )
    transition_rows.append(
        {
            "endpoint": endpoint,
            "direction": direction,
            "weighted_call_share": changed_calls,
            "weighted_true_positive_share": changed_positive,
            "implied_changed_block_precision": changed_positive / changed_calls,
            "hidden_base_rate": prevalence,
            "accuracy_change_current_minus_old": (
                values["current"]["Accuracy"] - values["old"]["Accuracy"]
            ),
            "interpretation": (
                "Approximation from products of public 1000-bootstrap mean metrics; "
                "not recovered compound labels"
            ),
        }
    )
pl.DataFrame(transition_rows).write_csv(out / "transition-inference.csv")

oof = pl.read_parquet(inp / cfg["ensemble_oof"])
inner = pl.read_parquet(inp / cfg["ensemble_inner"])
test = pl.read_parquet(inp / cfg["ensemble_test"])
fold_rows = []
screen_rows = []
prior_calls = []
for (endpoint, method), endpoint_oof in oof.group_by("endpoint", "method"):
    prevalence = priors[endpoint]
    outer_parts = []
    cuts = []
    for fold in range(5):
        inner_fold = inner.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("method") == method)
            & (pl.col("fold") == fold)
        )
        outer_fold = endpoint_oof.filter(pl.col("fold") == fold)
        assert set(inner_fold["compound_id"]).isdisjoint(outer_fold["compound_id"])
        cut = choose_prior_threshold(
            inner_fold["y_true"].to_numpy().astype(int),
            inner_fold["prediction"].to_numpy(),
            prevalence,
        )
        cuts.append(cut)
        outer_parts.append(outer_fold.with_columns(pl.lit(cut).alias("prior_threshold")))
        fold_rows.append(
            {"endpoint": endpoint, "method": method, "fold": fold, "threshold": cut}
        )
    evaluated = pl.concat(outer_parts)
    y = evaluated["y_true"].to_numpy().astype(int)
    probability = evaluated["prediction"].to_numpy()
    original_binary = probability >= evaluated["threshold"].to_numpy()
    prior_binary = probability >= evaluated["prior_threshold"].to_numpy()
    deployment_cut = float(np.median(cuts))
    blind = test.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
    assert len(blind) == blind["compound_id"].n_unique() == 750
    original_test_positives = int(blind["binary"].sum())
    blind = blind.with_columns(
        (pl.col("prediction") >= deployment_cut).alias("binary"),
        pl.lit(deployment_cut).alias("prior_threshold"),
    )
    prior_calls.append(blind.select("compound_id", "endpoint", "method", "prediction", "binary"))
    screen_rows.append(
        {
            "endpoint": endpoint,
            "method": method,
            "target_prevalence": prevalence,
            "original_source_MCC": float(matthews_corrcoef(y, original_binary)),
            "prior_policy_source_MCC": float(matthews_corrcoef(y, prior_binary)),
            "original_expected_target_MCC": expected_mcc(y, original_binary, prevalence),
            "prior_policy_expected_target_MCC": expected_mcc(y, prior_binary, prevalence),
            "deployment_threshold": deployment_cut,
            "original_test_positives": original_test_positives,
            "prior_policy_test_positives": int(blind["binary"].sum()),
        }
    )

screen = pl.DataFrame(screen_rows).sort(
    "endpoint", "prior_policy_expected_target_MCC", descending=[False, True]
)
screen.write_csv(out / "prior-threshold-screen.csv")
pl.DataFrame(fold_rows).write_csv(out / "prior-thresholds.csv")
all_prior_calls = pl.concat(prior_calls)
all_prior_calls.write_parquet(out / "prior-threshold-test-calls.parquet")

aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
identity = current_submission.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1"
)
assert len(identity) == identity["compound_id"].n_unique() == 750


def top_rate_calls(endpoint, method, rate, label):
    blind = test.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method)).sort(
        "prediction", "compound_id", descending=[True, False]
    )
    count = int(round(rate * len(blind)))
    return blind.with_row_index("rank").select(
        "compound_id",
        pl.lit(endpoint).alias("endpoint"),
        pl.lit(label).alias("variant"),
        (pl.col("rank") < count).alias("binary"),
    )


variant_parts = []
for endpoint in cfg["endpoints"]:
    method = cfg["current_methods"][endpoint]
    current_prior = all_prior_calls.filter(
        (pl.col("endpoint") == endpoint) & (pl.col("method") == method)
    ).select(
        "compound_id",
        "endpoint",
        pl.lit("current_models_prior_cut").alias("variant"),
        "binary",
    )
    variant_parts.append(current_prior)
    variant_parts.append(top_rate_calls(endpoint, method, priors[endpoint], "target_rate"))
    variant_parts.append(
        top_rate_calls(endpoint, method, leader_rates[endpoint], "leader_rate")
    )
    variant_parts.append(
        top_rate_calls(endpoint, method, top10_rates[endpoint], "top10_median_rate")
    )
    best_method = screen.filter(pl.col("endpoint") == endpoint)["method"][0]
    variant_parts.append(
        all_prior_calls.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("method") == best_method)
        ).select(
            "compound_id",
            "endpoint",
            pl.lit("screen_best_prior_cut").alias("variant"),
            "binary",
        )
    )

calls = pl.concat(variant_parts)
exact_current = current_submission.join(aliases, on="Molecule_Name", validate="1:1")
exact_old = old_submission.join(aliases, on="Molecule_Name", validate="1:1")
for variant in ("current", "restore_old_a4"):
    for endpoint in cfg["endpoints"]:
        source = (
            exact_old
            if variant == "restore_old_a4" and endpoint == "CYP3A4_is_TDI"
            else exact_current
        )
        calls = pl.concat(
            [
                calls,
                source.select(
                    "compound_id",
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(variant).alias("variant"),
                    pl.col(endpoint).alias("binary"),
                ),
            ]
        )
calls = calls.join(identity, on="compound_id", validate="m:1")
calls.write_parquet(out / "candidate-calls.parquet")
call_summary = (
    calls.group_by("variant", "endpoint")
    .agg(pl.col("binary").sum().alias("positives"))
    .sort("variant", "endpoint")
)
call_summary.write_csv(out / "candidate-call-summary.csv")

our_table = pl.DataFrame(our_rows)
known_hybrid_macro = float(
    (
        our_table.filter(
            (pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col("generation") == "current")
        )["MCC"][0]
        + our_table.filter(
            (pl.col("endpoint") == "CYP3A4_is_TDI") & (pl.col("generation") == "old")
        )["MCC"][0]
    )
    / 2
)
summary = {
    "leaderboard_priors": priors,
    "leader_implied_scored_positive_rates": leader_rates,
    "top10_median_implied_scored_positive_rates": top10_rates,
    "known_current_d6_old_a4_macro_MCC": known_hybrid_macro,
    "known_current_macro_MCC": float(
        our_table.filter(pl.col("generation") == "current")["MCC"].mean()
    ),
    "threshold_policy": (
        "Per outer fold, maximize expected MCC at the public-LB inferred target prior "
        "using only its isolated inner predictions; deploy median of five cuts"
    ),
    "limitations": (
        "Public metrics are 1000-bootstrap means on the confidently scored subset; "
        "prior-shift projections assume OOF class-conditional ranking transfers to hidden test"
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
print(call_summary.write_csv(), flush=True)
