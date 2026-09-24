"""Compare three submitted D6 decisions without inferring hidden labels."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key: str) -> pl.DataFrame:
    path = inputs / cfg[key]
    return pl.read_csv(path) if path.suffix == ".csv" else pl.read_parquet(path)


def package(key: str, alias: str) -> pl.DataFrame:
    table = read(key)
    if table.height != 750 or table["Molecule_Name"].n_unique() != 750:
        raise ValueError(f"Invalid 750-compound package: {key}")
    return table.select(
        pl.col("Molecule_Name").alias("molecule_name"),
        pl.col("SMILES").alias(f"{alias}_smiles"),
        pl.col(cfg["endpoint"]).cast(pl.Boolean).alias(f"{alias}_call"),
        pl.col("CYP3A4_is_TDI").cast(pl.Boolean).alias(f"{alias}_a4_call"),
    )


control = package("control_submission", "control")
support = package("support_submission", "support")
ap = package("ap_submission", "ap")
packages = control.join(support, on="molecule_name", validate="1:1").join(
    ap, on="molecule_name", validate="1:1"
)
if not (
    (packages["control_smiles"] == packages["support_smiles"]).all()
    and (packages["control_smiles"] == packages["ap_smiles"]).all()
    and (packages["control_a4_call"] == packages["support_a4_call"]).all()
    and (packages["control_a4_call"] == packages["ap_a4_call"]).all()
):
    raise ValueError("Package identities or CYP3A4 calls differ")
if [int(packages[f"{alias}_call"].sum()) for alias in ["control", "support", "ap"]] != [
    49,
    49,
    49,
]:
    raise ValueError("Unexpected submitted D6 positive counts")

test = read("test_diagnostics").join(
    packages.select("molecule_name", "control_call", "support_call", "ap_call"),
    on="molecule_name",
    validate="1:1",
)
if not (
    (test["old_call"] == test["control_call"]).all()
    and (test["new_call"] == test["support_call"]).all()
):
    raise ValueError("Historical support6 diagnostic does not match submitted vectors")
ap_test = read("ap_test").filter(pl.col("model") == cfg["ap_model"])
test = test.join(
    ap_test.select("compound_id", pl.col("call").alias("audited_ap_call")),
    on="compound_id",
    validate="1:1",
)
if not (test["ap_call"] == test["audited_ap_call"]).all():
    raise ValueError("AP-mass90 analytical and submitted calls differ")

test = test.with_columns(
    pl.concat_str(
        [
            pl.col(f"{alias}_call").cast(pl.Int8).cast(pl.String)
            for alias in ["control", "support", "ap"]
        ]
    ).alias("pattern")
)
test.write_parquet(outputs / "threeway-test.parquet")
numeric = [
    "top10_similarity_weighted_positive_rate",
    "nearest_positive_similarity",
    "nearest_negative_similarity",
    "nearest_any_similarity",
    "logp",
    "mw",
]
pattern_summary = (
    test.group_by("pattern")
    .agg(
        pl.len().alias("n"),
        *[pl.col(column).mean().alias(f"mean_{column}") for column in numeric],
    )
    .sort("pattern")
)
pattern_summary.write_csv(outputs / "blind-patterns.csv")

support_oof = read("support_oof")
ap_oof = read("ap_oof")


def oof_calls(source: pl.DataFrame, model: str, name: str) -> pl.DataFrame:
    return source.filter(pl.col("model") == model).select(
        "compound_id",
        pl.col("y_true").alias(f"{name}_label"),
        pl.col("call").alias(f"{name}_call"),
    )


oof = oof_calls(support_oof, cfg["control_model"], "control")
oof = oof.join(
    oof_calls(support_oof, cfg["support_model"], "support"),
    on="compound_id",
    validate="1:1",
).join(oof_calls(ap_oof, cfg["ap_model"], "ap"), on="compound_id", validate="1:1")
if oof.height != 1497 or not (
    (oof["control_label"] == oof["support_label"]).all()
    and (oof["control_label"] == oof["ap_label"]).all()
):
    raise ValueError("OOF keys or labels differ")
assignments = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
oof = oof.join(assignments, on="compound_id", validate="1:1").with_columns(
    pl.concat_str(
        [
            pl.col(f"{alias}_call").cast(pl.Int8).cast(pl.String)
            for alias in ["control", "support", "ap"]
        ]
    ).alias("pattern")
)
oof.write_parquet(outputs / "threeway-oof.parquet")
oof_pattern_summary = (
    oof.group_by("pattern")
    .agg(
        pl.len().alias("n"),
        pl.col("control_label").sum().alias("positives"),
        pl.col("control_label").mean().alias("positive_rate"),
    )
    .sort("pattern")
)
oof_pattern_summary.write_csv(outputs / "oof-patterns.csv")


def own_leaderboard(key: str) -> dict:
    row = read(key).filter(pl.col("Username") == cfg["username"])
    if row.height != 1:
        raise ValueError(f"Expected one team leaderboard row: {key}")
    return row.row(0, named=True)


leaderboards = {
    alias: own_leaderboard(f"{alias}_leaderboard") for alias in ["control", "support", "ap"]
}
leaderboard_rows = [
    {
        "model": alias,
        "submitted": row["Submitted"],
        "MCC": float(row["MCC"]),
        "precision": float(row["Precision"]),
        "recall": float(row["Recall"]),
        "F1": float(row["F1 Score"]),
        "estimated_TP_equivalents_per_750": float(row["Recall"])
        * cfg["estimated_hidden_prevalence"]
        * 750,
    }
    for alias, row in leaderboards.items()
]
pl.DataFrame(leaderboard_rows).write_csv(outputs / "leaderboard-comparison.csv")

rng = np.random.default_rng(cfg["permutation_seed"])


contrast_rows = []
for name, removed, added in [
    (
        "control_vs_support",
        test.filter(pl.col("control_call") & ~pl.col("support_call")),
        test.filter(~pl.col("control_call") & pl.col("support_call")),
    ),
    (
        "control_vs_ap",
        test.filter(pl.col("control_call") & ~pl.col("ap_call")),
        test.filter(~pl.col("control_call") & pl.col("ap_call")),
    ),
]:
    for feature in [
        "top10_similarity_weighted_positive_rate",
        "nearest_positive_similarity",
        "nearest_negative_similarity",
        "nearest_any_similarity",
    ]:
        a = removed[feature].to_numpy()
        b = added[feature].to_numpy()
        observed = float(a.mean() - b.mean())
        pooled = np.concatenate([a, b])
        extreme = 0
        for _ in range(cfg["permutation_replicates"]):
            shuffled = rng.permutation(pooled)
            delta = shuffled[: len(a)].mean() - shuffled[len(a) :].mean()
            extreme += abs(delta) >= abs(observed)
        contrast_rows.append(
            {
                "comparison": name,
                "feature": feature,
                "removed_n": len(a),
                "added_n": len(b),
                "removed_mean": float(a.mean()),
                "added_mean": float(b.mean()),
                "removed_minus_added": observed,
                "permutation_p_two_sided": (extreme + 1) / (cfg["permutation_replicates"] + 1),
            }
        )
pl.DataFrame(contrast_rows).write_csv(outputs / "blind-swap-contrasts.csv")

oof_swaps = []
for alias in ["support", "ap"]:
    removed = oof.filter(pl.col("control_call") & ~pl.col(f"{alias}_call"))
    added = oof.filter(~pl.col("control_call") & pl.col(f"{alias}_call"))
    if removed.height != added.height:
        raise ValueError("OOF call cohorts must be equal-sized")
    removed_positive = int(removed["control_label"].sum())
    added_positive = int(added["control_label"].sum())
    _, p = fisher_exact(
        [
            [added_positive, added.height - added_positive],
            [removed_positive, removed.height - removed_positive],
        ]
    )
    oof_swaps.append(
        {
            "candidate": alias,
            "swaps_each_direction": added.height,
            "added_positive": added_positive,
            "removed_positive": removed_positive,
            "net_positive": added_positive - removed_positive,
            "fisher_p_two_sided": float(p),
        }
    )
pl.DataFrame(oof_swaps).write_csv(outputs / "oof-swaps.csv")

summary = {
    "submitted_call_counts": {
        alias: int(packages[f"{alias}_call"].sum()) for alias in ["control", "support", "ap"]
    },
    "a4_exactly_unchanged": True,
    "test_pattern_counts": dict(zip(pattern_summary["pattern"], pattern_summary["n"], strict=True)),
    "oof_swap_counts": oof_swaps,
    "leaderboard_MCC": {row["model"]: row["MCC"] for row in leaderboard_rows},
    "limitations": [
        "Blind labels are unknown; aggregate LB metrics cannot identify individual errors.",
        "Reported MCC and recall are means over leaderboard bootstrap samples.",
        "Neighborhood comparisons are retrospective and do not validate a prospective selector.",
    ],
}
(outputs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
