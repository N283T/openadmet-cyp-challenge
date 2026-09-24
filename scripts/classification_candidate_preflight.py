"""Fail-closed preflight for fixed-call classification replacement candidates."""

import json
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import fisher_exact, spearmanr

from cyp_challenge.classification_preflight import (
    evaluate_preflight,
    exchangeable_effective_count,
)

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def decision(table: pl.DataFrame, model: str, prefix: str) -> pl.DataFrame:
    selected = table.filter(pl.col("model") == model)
    if selected.height == 0:
        raise ValueError(f"Missing decision model: {model}")
    columns = [
        "compound_id",
        pl.col("score").alias(f"{prefix}_score"),
        pl.col("call").alias(f"{prefix}_call"),
    ]
    if "y_true" in selected.columns:
        columns.insert(1, pl.col("y_true"))
    return selected.select(*columns).sort("compound_id")


def with_cohort(table: pl.DataFrame) -> pl.DataFrame:
    return table.with_columns(
        pl.when(pl.col("control_call") & pl.col("candidate_call"))
        .then(pl.lit("both_positive"))
        .when(pl.col("control_call"))
        .then(pl.lit("control_only"))
        .when(pl.col("candidate_call"))
        .then(pl.lit("candidate_only"))
        .otherwise(pl.lit("both_negative"))
        .alias("cohort")
    )


def pairwise_member_metrics(table: pl.DataFrame, domain: str) -> pl.DataFrame:
    members = cfg["member_models"]
    filtered = table.filter(pl.col("model").is_in(members))
    if set(filtered["model"].unique()) != set(members):
        raise ValueError(f"Incomplete {domain} member model set")
    value_column = "prediction" if "prediction" in filtered.columns else "score"
    scores = filtered.pivot(
        index="compound_id", on="model", values=value_column
    ).sort("compound_id")
    calls = filtered.pivot(index="compound_id", on="model", values="call").sort(
        "compound_id"
    )
    rows = []
    for left, right in combinations(members, 2):
        left_call = calls[left].to_numpy().astype(bool)
        right_call = calls[right].to_numpy().astype(bool)
        union = int((left_call | right_call).sum())
        rows.append(
            {
                "domain": domain,
                "model_a": left,
                "model_b": right,
                "spearman": float(spearmanr(scores[left], scores[right]).statistic),
                "call_phi": float(
                    np.corrcoef(left_call.astype(float), right_call.astype(float))[0, 1]
                ),
                "call_jaccard": float((left_call & right_call).sum() / union),
            }
        )
    return pl.DataFrame(rows)


decision_test = read("decision_test")
test = decision(decision_test, cfg["control_model"], "control").join(
    decision(decision_test, cfg["candidate_model"], "candidate"),
    on="compound_id",
    validate="1:1",
)
test = with_cohort(test)
control_calls = int(test["control_call"].sum())
candidate_calls = int(test["candidate_call"].sum())
expected_calls = int(cfg["expected_positive_calls"])
if control_calls != expected_calls or candidate_calls != expected_calls:
    raise ValueError(
        f"Expected {expected_calls} calls; control={control_calls}, candidate={candidate_calls}"
    )
control_only = test.filter(pl.col("cohort") == "control_only")
candidate_only = test.filter(pl.col("cohort") == "candidate_only")
if control_only.height != candidate_only.height or control_only.height == 0:
    raise ValueError("Preflight requires equal, nonzero fixed-call swap cohorts")

decision_oof = read("decision_oof")
control_oof = decision(decision_oof, cfg["control_model"], "control")
candidate_oof = decision(decision_oof, cfg["candidate_model"], "candidate").drop(
    "y_true"
)
oof = with_cohort(control_oof.join(candidate_oof, on="compound_id", validate="1:1"))
oof_summary = oof.group_by("cohort").agg(
    pl.len().alias("n"),
    pl.col("y_true").sum().alias("positives"),
    pl.col("y_true").mean().alias("positive_rate"),
).sort("cohort")
oof_summary.write_csv(out / "oof-swap-summary.csv")
oof_by_cohort = {row["cohort"]: row for row in oof_summary.iter_rows(named=True)}
oof_control = oof_by_cohort["control_only"]
oof_candidate = oof_by_cohort["candidate_only"]
_, oof_p = fisher_exact(
    [
        [
            oof_candidate["positives"],
            oof_candidate["n"] - oof_candidate["positives"],
        ],
        [oof_control["positives"], oof_control["n"] - oof_control["positives"]],
    ]
)

dependence = pl.concat(
    [
        pairwise_member_metrics(read("member_oof"), "oof"),
        pairwise_member_metrics(read("member_test"), "test"),
    ]
)
dependence.write_csv(out / "member-dependence.csv")
test_dependence = dependence.filter(pl.col("domain") == "test")
oof_dependence = dependence.filter(pl.col("domain") == "oof")
mean_test_spearman = float(test_dependence["spearman"].mean())
mean_oof_spearman = float(oof_dependence["spearman"].mean())
mean_test_phi = float(test_dependence["call_phi"].mean())
effective_votes = exchangeable_effective_count(len(cfg["member_models"]), mean_test_phi)

generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=int(cfg["morgan_radius"]), fpSize=int(cfg["morgan_bits"])
)
compounds = read("compounds")
labels = (
    read("measurements")
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", pl.col("value").cast(pl.Int8).alias("y_true"))
    .join(compounds, on="compound_id", validate="1:1")
)
train_y = labels["y_true"].to_numpy()
train_fps = []
for smiles in labels["standardized_smiles"]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid training SMILES: {smiles}")
    train_fps.append(generator.GetFingerprint(molecule))

swapped = test.filter(pl.col("cohort").is_in(["control_only", "candidate_only"])).join(
    compounds, on="compound_id", validate="1:1"
)
neighbor_rows = []
for row in swapped.iter_rows(named=True):
    molecule = Chem.MolFromSmiles(row["standardized_smiles"])
    if molecule is None:
        raise ValueError(f"Invalid test SMILES: {row['standardized_smiles']}")
    fingerprint = generator.GetFingerprint(molecule)
    similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fingerprint, train_fps))
    order = np.argsort(-similarities, kind="stable")[: int(cfg["knn_k"])]
    weights = similarities[order]
    neighbor_rows.append(
        {
            "compound_id": row["compound_id"],
            "cohort": row["cohort"],
            "nearest_any_similarity": float(similarities.max()),
            "nearest_positive_similarity": float(similarities[train_y == 1].max()),
            "nearest_negative_similarity": float(similarities[train_y == 0].max()),
            "top10_weighted_positive_rate": float(
                np.average(train_y[order], weights=weights)
            ),
        }
    )
neighbors = pl.DataFrame(neighbor_rows)
neighbors.write_csv(out / "swapped-neighborhoods.csv")

control_neighbor = neighbors.filter(pl.col("cohort") == "control_only")[
    "top10_weighted_positive_rate"
].to_numpy()
candidate_neighbor = neighbors.filter(pl.col("cohort") == "candidate_only")[
    "top10_weighted_positive_rate"
].to_numpy()
neighbor_delta = float(control_neighbor.mean() - candidate_neighbor.mean())
combined = np.concatenate([control_neighbor, candidate_neighbor])
rng = np.random.default_rng(int(cfg["random_seed"]))
permutation_delta = np.empty(int(cfg["permutation_replicates"]))
for index in range(len(permutation_delta)):
    shuffled = rng.permutation(combined)
    permutation_delta[index] = shuffled[: len(control_neighbor)].mean() - shuffled[
        len(control_neighbor) :
    ].mean()
neighbor_p = float(
    (np.count_nonzero(np.abs(permutation_delta) >= abs(neighbor_delta)) + 1)
    / (len(permutation_delta) + 1)
)

candidate_vs_control_spearman = float(
    spearmanr(test["candidate_score"], test["control_score"]).statistic
)
metrics = {
    "positive_calls": expected_calls,
    "swaps_each_direction": control_only.height,
    "swap_fraction": control_only.height / expected_calls,
    "candidate_vs_control_test_spearman": candidate_vs_control_spearman,
    "mean_test_member_score_spearman": mean_test_spearman,
    "mean_oof_member_score_spearman": mean_oof_spearman,
    "test_minus_oof_member_spearman": mean_test_spearman - mean_oof_spearman,
    "mean_test_member_call_phi": mean_test_phi,
    "effective_vote_count": effective_votes,
    "oof_candidate_minus_control_swap_positive_rate": float(
        oof_candidate["positive_rate"] - oof_control["positive_rate"]
    ),
    "oof_swap_fisher_p_two_sided": float(oof_p),
    "control_minus_candidate_top10_weighted_positive_rate": neighbor_delta,
    "top10_weighted_positive_rate_permutation_p": neighbor_p,
}
status, checks = evaluate_preflight(metrics, cfg["thresholds"])
pl.DataFrame(checks).write_csv(out / "checks.csv")
summary = {
    "status": status,
    "endpoint": cfg["endpoint"],
    "control_model": cfg["control_model"],
    "candidate_model": cfg["candidate_model"],
    "metrics": metrics,
    "thresholds": cfg["thresholds"],
    "checks": checks,
    "enforced": bool(cfg.get("enforce", True)),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
if summary["enforced"] and status != "PASS":
    raise SystemExit(2)
