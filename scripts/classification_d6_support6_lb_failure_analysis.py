"""Diagnose why the support6 CYP2D6 consensus failed on the leaderboard."""

import json
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator
from scipy.stats import fisher_exact, spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def parquet(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def csv(key: str) -> pl.DataFrame:
    return pl.read_csv(inp / cfg[key])


def own_row(table: pl.DataFrame) -> dict:
    selected = table.filter(pl.col("Username") == cfg["our_username"])
    assert selected.height == 1
    return selected.row(0, named=True)


def selected(table: pl.DataFrame, model: str, prefix: str) -> pl.DataFrame:
    return table.filter(pl.col("model") == model).select(
        "compound_id",
        pl.col("score").alias(f"{prefix}_score"),
        pl.col("call").alias(f"{prefix}_call"),
    )


pre = own_row(csv("pre_d6_leaderboard"))
post = own_row(csv("post_d6_leaderboard"))
prevalence = cfg["hidden_prevalence"]
leaderboard_rows = []
for generation, row in [("control", pre), ("support6", post)]:
    call_rate = float(row["Recall"]) * prevalence / float(row["Precision"])
    tp_share = float(row["Recall"]) * prevalence
    leaderboard_rows.append(
        {
            "generation": generation,
            "rank": int(row["Rank"]),
            "submitted": row["Submitted"],
            "MCC": float(row["MCC"]),
            "accuracy": float(row["Accuracy"]),
            "precision": float(row["Precision"]),
            "recall": float(row["Recall"]),
            "F1": float(row["F1 Score"]),
            "implied_call_rate": call_rate,
            "estimated_tp_share": tp_share,
            "estimated_fp_share": call_rate - tp_share,
        }
    )
leaderboard = pl.DataFrame(leaderboard_rows)
leaderboard.write_csv(out / "leaderboard-movement.csv")

control = leaderboard.filter(pl.col("generation") == "control").row(0, named=True)
candidate = leaderboard.filter(pl.col("generation") == "support6").row(0, named=True)
lb_delta = {
    field: candidate[field] - control[field]
    for field in [
        "MCC",
        "accuracy",
        "precision",
        "recall",
        "F1",
        "implied_call_rate",
        "estimated_tp_share",
        "estimated_fp_share",
    ]
}
lb_delta["rank"] = candidate["rank"] - control["rank"]
lb_delta["estimated_net_true_positives_lost_per_750"] = -750 * lb_delta[
    "estimated_tp_share"
]

consensus_test = parquet("consensus_test")
test = selected(consensus_test, cfg["deployed_model"], "old").join(
    selected(consensus_test, cfg["candidate_model"], "new"),
    on="compound_id",
    validate="1:1",
)
test = test.with_columns(
    pl.when(pl.col("old_call") & pl.col("new_call"))
    .then(pl.lit("both_positive"))
    .when(pl.col("old_call"))
    .then(pl.lit("old_only"))
    .when(pl.col("new_call"))
    .then(pl.lit("new_only"))
    .otherwise(pl.lit("both_negative"))
    .alias("cohort")
)

members = [
    f"{root}::{split}"
    for root in cfg["member_roots"]
    for split in cfg["member_splits"]
]
member_test = parquet("member_test").filter(pl.col("model").is_in(members))
member_votes = member_test.group_by("compound_id").agg(
    pl.col("call").sum().alias("support6_votes"),
    pl.col("prediction").mean().alias("support6_raw_mean"),
    pl.col("call")
    .filter(pl.col("model").str.ends_with("::random"))
    .sum()
    .alias("random_votes"),
    pl.col("call")
    .filter(pl.col("model").str.ends_with("::scaffold"))
    .sum()
    .alias("scaffold_votes"),
)
test = test.join(member_votes, on="compound_id", validate="1:1")

for model in ["support6_equal_rank", "all9_top_call_vote", "all9_equal_rank"]:
    extra = consensus_test.filter(pl.col("model") == model).select(
        "compound_id", pl.col("call").alias(f"{model}_call")
    )
    test = test.join(extra, on="compound_id", validate="1:1")

test = (
    test.join(parquet("test_support"), on="compound_id", validate="1:1")
    .join(
        parquet("test_mapping").select("compound_id", "molecule_name"),
        on="compound_id",
        validate="1:1",
    )
    .join(parquet("compounds"), on="compound_id", validate="1:1")
)

generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg["morgan_radius"], fpSize=cfg["morgan_bits"]
)


def molecule_features(smiles: str) -> tuple:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return (
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        Lipinski.RingCount(mol),
        Lipinski.FractionCSP3(mol),
        generator.GetFingerprint(mol),
    )


measurements = parquet("measurements").filter(pl.col("endpoint") == cfg["endpoint"])
train = measurements.select("compound_id", pl.col("value").cast(pl.Int8).alias("y_true")).join(
    parquet("compounds"), on="compound_id", validate="1:1"
)
train_features = [molecule_features(s) for s in train["standardized_smiles"]]
train_fps = [row[-1] for row in train_features]
train_y = train["y_true"].to_numpy()

feature_rows = []
for row in test.select("compound_id", "standardized_smiles").iter_rows(named=True):
    mw, logp, tpsa, hbd, hba, rings, fraction_csp3, fp = molecule_features(
        row["standardized_smiles"]
    )
    similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, train_fps))
    order = np.argsort(-similarities, kind="stable")
    top = order[: cfg["knn_k"]]
    top5 = order[:5]
    positive_similarities = similarities[train_y == 1]
    negative_similarities = similarities[train_y == 0]
    weights = similarities[top]
    feature_rows.append(
        {
            "compound_id": row["compound_id"],
            "mw": mw,
            "logp": logp,
            "tpsa": tpsa,
            "hbd": hbd,
            "hba": hba,
            "rings": rings,
            "fraction_csp3": fraction_csp3,
            "nearest_any_similarity": float(similarities[order[0]]),
            "nearest_positive_similarity": float(positive_similarities.max()),
            "nearest_negative_similarity": float(negative_similarities.max()),
            "nearest_neighbor_label": int(train_y[order[0]]),
            "top5_positive_rate": float(train_y[top5].mean()),
            "top10_positive_rate": float(train_y[top].mean()),
            "top10_similarity_weighted_positive_rate": float(
                np.average(train_y[top], weights=weights)
            ),
        }
    )
test = test.join(pl.DataFrame(feature_rows), on="compound_id", validate="1:1")
assert float((test["support"] - test["nearest_any_similarity"]).abs().max()) < 1e-12

numeric = [
    "old_score",
    "new_score",
    "support6_votes",
    "random_votes",
    "scaffold_votes",
    "support",
    "mw",
    "logp",
    "tpsa",
    "rings",
    "fraction_csp3",
    "nearest_any_similarity",
    "nearest_positive_similarity",
    "nearest_negative_similarity",
    "nearest_neighbor_label",
    "top5_positive_rate",
    "top10_positive_rate",
    "top10_similarity_weighted_positive_rate",
]
aggregations = [pl.len().alias("n")]
for column in numeric:
    aggregations.extend(
        [
            pl.col(column).mean().alias(f"mean_{column}"),
            pl.col(column).median().alias(f"median_{column}"),
        ]
    )
aggregations.extend(
    [
        pl.col("support6_equal_rank_call").sum().alias("support6_equal_rank_calls"),
        pl.col("all9_top_call_vote_call").sum().alias("all9_vote_calls"),
        pl.col("all9_equal_rank_call").sum().alias("all9_equal_rank_calls"),
    ]
)
test_summary = test.group_by("cohort").agg(aggregations).sort("cohort")
test_summary.write_csv(out / "test-swap-cohort-summary.csv")
test.filter(pl.col("cohort").is_in(["old_only", "new_only"])).sort(
    "cohort", "support6_votes", descending=[False, True]
).write_csv(out / "swapped-compounds.csv")
test.write_parquet(out / "test-diagnostics.parquet")

rng = np.random.default_rng(42)
contrast_rows = []
old_swap = test.filter(pl.col("cohort") == "old_only")
new_swap = test.filter(pl.col("cohort") == "new_only")
for column in [
    "support",
    "nearest_positive_similarity",
    "nearest_negative_similarity",
    "nearest_neighbor_label",
    "top5_positive_rate",
    "top10_positive_rate",
    "top10_similarity_weighted_positive_rate",
    "logp",
]:
    old_values = old_swap[column].to_numpy().astype(float)
    new_values = new_swap[column].to_numpy().astype(float)
    observed = float(old_values.mean() - new_values.mean())
    combined = np.concatenate([old_values, new_values])
    permutation_deltas = np.empty(10000)
    for index in range(len(permutation_deltas)):
        shuffled = rng.permutation(combined)
        permutation_deltas[index] = shuffled[: len(old_values)].mean() - shuffled[
            len(old_values) :
        ].mean()
    bootstrap_deltas = np.empty(10000)
    for index in range(len(bootstrap_deltas)):
        bootstrap_deltas[index] = rng.choice(old_values, len(old_values)).mean() - rng.choice(
            new_values, len(new_values)
        ).mean()
    contrast_rows.append(
        {
            "metric": column,
            "old_only_mean": float(old_values.mean()),
            "new_only_mean": float(new_values.mean()),
            "old_minus_new": observed,
            "bootstrap_low": float(np.quantile(bootstrap_deltas, 0.025)),
            "bootstrap_high": float(np.quantile(bootstrap_deltas, 0.975)),
            "permutation_p_two_sided": float(
                (np.count_nonzero(np.abs(permutation_deltas) >= abs(observed)) + 1)
                / (len(permutation_deltas) + 1)
            ),
        }
    )
pl.DataFrame(contrast_rows).write_csv(out / "test-swap-contrasts.csv")

consensus_oof = parquet("consensus_oof")
old_oof = consensus_oof.filter(pl.col("model") == cfg["deployed_model"]).select(
    "compound_id",
    "y_true",
    pl.col("score").alias("old_score"),
    pl.col("call").alias("old_call"),
)
new_oof = consensus_oof.filter(pl.col("model") == cfg["candidate_model"]).select(
    "compound_id", pl.col("score").alias("new_score"), pl.col("call").alias("new_call")
)
oof = old_oof.join(new_oof, on="compound_id", validate="1:1").join(
    parquet("oof_support").select("compound_id", "support", "support_bin"),
    on="compound_id",
    validate="1:1",
)
oof = oof.with_columns(
    pl.when(pl.col("old_call") & pl.col("new_call"))
    .then(pl.lit("both_positive"))
    .when(pl.col("old_call"))
    .then(pl.lit("old_only"))
    .when(pl.col("new_call"))
    .then(pl.lit("new_only"))
    .otherwise(pl.lit("both_negative"))
    .alias("cohort")
)
oof_summary = oof.group_by("cohort").agg(
    pl.len().alias("n"),
    pl.col("y_true").sum().alias("positives"),
    pl.col("y_true").mean().alias("positive_rate"),
    pl.col("support").mean().alias("mean_support"),
    pl.col("support").median().alias("median_support"),
    pl.col("old_score").mean().alias("mean_old_score"),
    pl.col("new_score").mean().alias("mean_new_score"),
).sort("cohort")
oof_summary.write_csv(out / "oof-swap-cohort-summary.csv")
oof.group_by("cohort", "support_bin").agg(
    pl.len().alias("n"),
    pl.col("y_true").sum().alias("positives"),
    pl.col("y_true").mean().alias("positive_rate"),
).sort("cohort", "support_bin").write_csv(out / "oof-swap-by-support-bin.csv")


def dependence(table: pl.DataFrame, domain: str) -> list[dict]:
    filtered = table.filter(pl.col("model").is_in(members)).select(
        "compound_id", "model", "call"
    )
    wide = filtered.pivot(index="compound_id", on="model", values="call").sort("compound_id")
    rows = []
    for left, right in combinations(members, 2):
        a = wide[left].to_numpy().astype(bool)
        b = wide[right].to_numpy().astype(bool)
        union = int((a | b).sum())
        rows.append(
            {
                "domain": domain,
                "model_a": left,
                "model_b": right,
                "phi": float(np.corrcoef(a.astype(float), b.astype(float))[0, 1]),
                "jaccard": float((a & b).sum() / union) if union else 1.0,
                "overlap": int((a & b).sum()),
                "union": union,
            }
        )
    return rows


dependence_rows = dependence(parquet("member_test"), "test") + dependence(
    parquet("member_oof"), "oof"
)
dependence_frame = pl.DataFrame(dependence_rows)
dependence_frame.write_csv(out / "member-call-dependence.csv")


def score_dependence(table: pl.DataFrame, domain: str) -> list[dict]:
    filtered = table.filter(
        pl.col("model").is_in([*members, cfg["deployed_model"]])
    ).select("compound_id", "model", "prediction")
    wide = filtered.pivot(index="compound_id", on="model", values="prediction").sort(
        "compound_id"
    )
    rows = []
    for left, right in combinations(members, 2):
        rows.append(
            {
                "domain": domain,
                "comparison": "within_support6",
                "model_a": left,
                "model_b": right,
                "spearman": float(spearmanr(wide[left], wide[right]).statistic),
            }
        )
    for member in members:
        rows.append(
            {
                "domain": domain,
                "comparison": "member_vs_deployed",
                "model_a": member,
                "model_b": cfg["deployed_model"],
                "spearman": float(
                    spearmanr(wide[member], wide[cfg["deployed_model"]]).statistic
                ),
            }
        )
    return rows


score_dependence_frame = pl.DataFrame(
    score_dependence(parquet("member_test"), "test")
    + score_dependence(parquet("member_oof"), "oof")
)
score_dependence_frame.write_csv(out / "member-score-dependence.csv")

composite_wide = consensus_test.pivot(
    index="compound_id", on="model", values="score"
).sort("compound_id")
composite_calls = consensus_test.pivot(
    index="compound_id", on="model", values="call"
).sort("compound_id")
composite_rows = []
for model in [
    "support6_top_call_vote",
    "support6_equal_rank",
    "all9_top_call_vote",
    "all9_equal_rank",
]:
    left = composite_calls[model].to_numpy().astype(bool)
    right = composite_calls[cfg["deployed_model"]].to_numpy().astype(bool)
    composite_rows.append(
        {
            "model": model,
            "spearman_vs_deployed": float(
                spearmanr(
                    composite_wide[model], composite_wide[cfg["deployed_model"]]
                ).statistic
            ),
            "call_phi_vs_deployed": float(
                np.corrcoef(left.astype(float), right.astype(float))[0, 1]
            ),
            "call_jaccard_vs_deployed": float((left & right).sum() / (left | right).sum()),
            "positive_call_overlap": int((left & right).sum()),
        }
    )
pl.DataFrame(composite_rows).write_csv(out / "composite-dependence.csv")

swap_counts = test.group_by("cohort").len().to_dict(as_series=False)
swap_counts = dict(zip(swap_counts["cohort"], swap_counts["len"], strict=True))
oof_rates = {
    row["cohort"]: row["positive_rate"] for row in oof_summary.iter_rows(named=True)
}
test_dependence = dependence_frame.filter(pl.col("domain") == "test")
oof_new = oof_summary.filter(pl.col("cohort") == "new_only").row(0, named=True)
oof_old = oof_summary.filter(pl.col("cohort") == "old_only").row(0, named=True)
_, oof_fisher_p = fisher_exact(
    [
        [oof_new["positives"], oof_new["n"] - oof_new["positives"]],
        [oof_old["positives"], oof_old["n"] - oof_old["positives"]],
    ]
)
mean_test_phi = float(test_dependence["phi"].mean())
test_score_dependence = score_dependence_frame.filter(
    (pl.col("domain") == "test") & (pl.col("comparison") == "within_support6")
)
oof_score_dependence = score_dependence_frame.filter(
    (pl.col("domain") == "oof") & (pl.col("comparison") == "within_support6")
)
summary = {
    "leaderboard_delta": lb_delta,
    "test_call_transitions": swap_counts,
    "oof_swap_positive_rates": oof_rates,
    "oof_new_minus_old_only_positive_rate": oof_rates["new_only"]
    - oof_rates["old_only"],
    "oof_new_minus_old_only_positive_count": int(
        oof_new["positives"] - oof_old["positives"]
    ),
    "oof_swap_fisher_exact_p_two_sided": float(oof_fisher_p),
    "mean_test_member_phi": mean_test_phi,
    "mean_test_member_jaccard": float(test_dependence["jaccard"].mean()),
    "mean_test_member_score_spearman": float(test_score_dependence["spearman"].mean()),
    "mean_oof_member_score_spearman": float(oof_score_dependence["spearman"].mean()),
    "exchangeable_correlation_effective_vote_count": 6
    / (1 + 5 * mean_test_phi),
    "limitations": [
        "Leaderboard metrics are bootstrap means on confidently scored hidden labels; "
        "implied shares are weighted estimates, not recovered molecule labels.",
        "Chemical-neighbor diagnostics are descriptive and cannot identify which "
        "individual swapped compounds are hidden positives.",
    ],
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
