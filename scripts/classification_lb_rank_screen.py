"""Evaluate saved classifiers at public-LB-derived positive-call rates."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
lb = json.loads((inp / cfg["lb_summary"]).read_text())
our_lb = pl.read_csv(inp / cfg["our_lb"])
leaderboard = pl.read_csv(inp / cfg["leaderboard_points"])


def expected_mcc(y, binary, prevalence):
    positive = y == 1
    tpr = float(binary[positive].mean())
    fpr = float(binary[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator > 0 else 0.0


def rank_scores(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks / max(len(values) - 1, 1)


def top_rate(values, rate):
    count = int(round(rate * len(values)))
    order = np.argsort(-values, kind="stable")
    result = np.zeros(len(values), dtype=bool)
    result[order[:count]] = True
    return result


ensemble_oof = pl.read_parquet(inp / cfg["ensemble_oof"])
ensemble_test = pl.read_parquet(inp / cfg["ensemble_test"])
caruana_oof = pl.read_parquet(inp / cfg["caruana_oof"])
caruana_test = pl.read_parquet(inp / cfg["caruana_test"])
weighted_oof = pl.read_parquet(inp / cfg["weighted_oof"])
weighted_test = pl.read_parquet(inp / cfg["weighted_test"])


def component(source, method, endpoint, test=False):
    table = {
        ("ensemble", False): ensemble_oof,
        ("ensemble", True): ensemble_test,
        ("caruana", False): caruana_oof,
        ("caruana", True): caruana_test,
        ("weighted", False): weighted_oof,
        ("weighted", True): weighted_test,
    }[(source, test)]
    result = table.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
    if test and result["compound_id"].n_unique() != len(result):
        result = result.group_by("compound_id", "endpoint", "method").agg(
            pl.col("prediction").mean()
        )
    expected = 750 if test else (1497 if endpoint == "CYP2D6_is_TDI" else 3584)
    assert len(result) == result["compound_id"].n_unique() == expected
    return result.sort("compound_id")


definitions = {
    "CYP2D6_is_TDI": {
        "current": [("ensemble", "all_plain_pfn_graph_mean", 1.0)],
        "balanced_lgbm": [("ensemble", "lgbm_balanced_mean", 1.0)],
        "caruana_auc": [("caruana", "caruana_auc", 1.0)],
        "current_caruana_rank_mean": [
            ("ensemble", "all_plain_pfn_graph_mean", 0.5),
            ("caruana", "caruana_auc", 0.5),
        ],
        "current_balanced_rank_mean": [
            ("ensemble", "all_plain_pfn_graph_mean", 0.5),
            ("ensemble", "lgbm_balanced_mean", 0.5),
        ],
        "lb_weighted_mean": [("weighted", "six_view_mean", 1.0)],
        "current_lb_weighted_rank_mean": [
            ("ensemble", "all_plain_pfn_graph_mean", 0.5),
            ("weighted", "six_view_mean", 0.5),
        ],
        "current_caruana_lb_weighted_rank_mean": [
            ("ensemble", "all_plain_pfn_graph_mean", 1 / 3),
            ("caruana", "caruana_auc", 1 / 3),
            ("weighted", "six_view_mean", 1 / 3),
        ],
    },
    "CYP3A4_is_TDI": {
        "current": [("ensemble", "lgbm_balanced_mean", 1.0)],
        "historical": [("caruana", "fixed_historical", 1.0)],
        "caruana_auc": [("caruana", "caruana_auc", 1.0)],
        "historical_caruana_rank_mean": [
            ("caruana", "fixed_historical", 0.5),
            ("caruana", "caruana_auc", 0.5),
        ],
        "historical_current_rank_mean": [
            ("caruana", "fixed_historical", 0.5),
            ("ensemble", "lgbm_balanced_mean", 0.5),
        ],
    },
}


def build_candidate(endpoint, members, test=False):
    parts = [component(source, method, endpoint, test=test) for source, method, _ in members]
    keys = ["compound_id"] + ([] if test else ["fold", "y_true"])
    reference = parts[0].select(keys)
    for part in parts[1:]:
        assert part.select(keys).equals(reference)
    if len(parts) == 1:
        score = parts[0]["prediction"].to_numpy()
    else:
        score = np.zeros(len(parts[0]), dtype=float)
        if test:
            for part, (_, _, weight) in zip(parts, members, strict=True):
                score += weight * rank_scores(part["prediction"].to_numpy())
        else:
            folds = parts[0]["fold"].to_numpy()
            for fold in range(5):
                mask = folds == fold
                for part, (_, _, weight) in zip(parts, members, strict=True):
                    score[mask] += weight * rank_scores(part["prediction"].to_numpy()[mask])
    return reference.with_columns(pl.Series("score", score))


rates = {}
for endpoint in cfg["endpoints"]:
    old_rate = float(
        our_lb.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("generation") == "old")
        )["implied_scored_positive_rate"][0]
    )
    current_rate = float(
        our_lb.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("generation") == "current")
        )["implied_scored_positive_rate"][0]
    )
    rates[endpoint] = {
        "target_prior": lb["leaderboard_priors"][endpoint],
        "leader_rate": lb["leader_implied_scored_positive_rates"][endpoint],
        "top10_median_rate": lb["top10_median_implied_scored_positive_rates"][endpoint],
        "our_old_scored_rate": old_rate,
        "our_current_scored_rate": current_rate,
        "old_full_rate" if endpoint == "CYP3A4_is_TDI" else "current_full_rate": (
            249 / 750 if endpoint == "CYP3A4_is_TDI" else 43 / 750
        ),
    }

screen_rows = []
test_call_parts = []
for endpoint, candidates in definitions.items():
    prevalence = lb["leaderboard_priors"][endpoint]
    for candidate_name, members in candidates.items():
        candidate_oof = build_candidate(endpoint, members)
        candidate_test = build_candidate(endpoint, members, test=True)
        y = candidate_oof["y_true"].to_numpy().astype(int)
        folds = candidate_oof["fold"].to_numpy()
        score = candidate_oof["score"].to_numpy()
        for policy, rate in rates[endpoint].items():
            binary = np.zeros(len(y), dtype=bool)
            for fold in range(5):
                mask = folds == fold
                binary[mask] = top_rate(score[mask], rate)
            blind_binary = top_rate(candidate_test["score"].to_numpy(), rate)
            screen_rows.append(
                {
                    "endpoint": endpoint,
                    "candidate": candidate_name,
                    "rate_policy": policy,
                    "rate": rate,
                    "source_MCC": float(matthews_corrcoef(y, binary)),
                    "source_AUROC": float(roc_auc_score(y, score)),
                    "expected_target_MCC": expected_mcc(y, binary, prevalence),
                    "test_positives": int(blind_binary.sum()),
                }
            )
            test_call_parts.append(
                candidate_test.select("compound_id").with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(candidate_name).alias("candidate"),
                    pl.lit(policy).alias("rate_policy"),
                    pl.Series("binary", blind_binary),
                )
            )

screen = pl.DataFrame(screen_rows).sort(
    "endpoint", "rate_policy", "expected_target_MCC", descending=[False, False, True]
)
screen.write_csv(out / "rank-rate-screen.csv")
test_calls = pl.concat(test_call_parts)
test_calls.write_parquet(out / "rank-rate-test-calls.parquet")

best_rows = (
    screen.sort(
        "endpoint", "rate_policy", "expected_target_MCC", descending=[False, False, True]
    )
    .group_by("endpoint", "rate_policy", maintain_order=True)
    .first()
)
best_rows.write_csv(out / "best-by-rate.csv")

old_submission = pl.read_parquet(inp / cfg["old_submission"]).sort("Molecule_Name")
current_submission = pl.read_parquet(inp / cfg["current_submission"]).sort("Molecule_Name")
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
identity = current_submission.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1"
)

specific = {
    "best_top10_rate": {
        endpoint: (
            best_rows.filter(
                (pl.col("endpoint") == endpoint)
                & (pl.col("rate_policy") == "top10_median_rate")
            )["candidate"][0],
            "top10_median_rate",
        )
        for endpoint in cfg["endpoints"]
    },
    "current_d6_caruana_a4_old_count": {
        "CYP2D6_is_TDI": ("exact_current", "exact"),
        "CYP3A4_is_TDI": ("caruana_auc", "old_full_rate"),
    },
    "current_d6_rankblend_a4_old_count": {
        "CYP2D6_is_TDI": ("exact_current", "exact"),
        "CYP3A4_is_TDI": ("historical_caruana_rank_mean", "old_full_rate"),
    },
}
wide_rows = []
comparison_rows = []
for variant, endpoint_choices in specific.items():
    parts = []
    for endpoint, (candidate_name, policy) in endpoint_choices.items():
        if policy == "exact":
            values = current_submission.join(aliases, on="Molecule_Name", validate="1:1").select(
                "compound_id", pl.col(endpoint).alias(endpoint)
            )
        else:
            values = test_calls.filter(
                (pl.col("endpoint") == endpoint)
                & (pl.col("candidate") == candidate_name)
                & (pl.col("rate_policy") == policy)
            ).select("compound_id", pl.col("binary").alias(endpoint))
        parts.append(values)
    wide = identity.join(parts[0], on="compound_id", validate="1:1").join(
        parts[1], on="compound_id", validate="1:1"
    )
    wide = wide.select("SMILES", "Molecule_Name", *cfg["endpoints"]).sort("Molecule_Name")
    wide.write_parquet(out / f"{variant}.parquet")
    validation = validate(
        out / f"{variant}.parquet", "classification", inp / cfg["current_submission"]
    )
    wide_rows.append(
        {
            "variant": variant,
            **{endpoint: int(wide[endpoint].sum()) for endpoint in cfg["endpoints"]},
            "valid": validation["valid"],
        }
    )
    for endpoint in cfg["endpoints"]:
        reference = current_submission if endpoint == "CYP2D6_is_TDI" else old_submission
        reference_binary = reference[endpoint].to_numpy().astype(bool)
        variant_binary = wide[endpoint].to_numpy().astype(bool)
        comparison_rows.append(
            {
                "variant": variant,
                "endpoint": endpoint,
                "reference": "current_D6" if endpoint == "CYP2D6_is_TDI" else "old_A4",
                "both_positive": int((reference_binary & variant_binary).sum()),
                "reference_only": int((reference_binary & ~variant_binary).sum()),
                "variant_only": int((~reference_binary & variant_binary).sum()),
                "both_negative": int((~reference_binary & ~variant_binary).sum()),
            }
        )
pl.DataFrame(comparison_rows).write_csv(out / "variant-reference-overlap.csv")

leader_rows = {
    endpoint: leaderboard.filter(
        (pl.col("endpoint") == endpoint) & (pl.col("rank") == 1)
    ).row(0, named=True)
    for endpoint in cfg["endpoints"]
}
summary = {
    "rates": rates,
    "variants": wide_rows,
    "leader_rows": leader_rows,
    "evaluation": (
        "Fixed positive-call rates are external LB facts; within each outer fold, calls use "
        "prediction rank only. Expected target MCC substitutes the inferred hidden prior "
        "while retaining OOF TPR/FPR. Candidate choice across screened models remains exploratory."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
