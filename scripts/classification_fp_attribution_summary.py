"""Summarize robust, interpretable rows from the Morgan attribution table."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
bits = pl.read_csv(inp / cfg["attributions"])
ablations = pl.read_csv(inp / cfg["ablations"])

selected_tables = []
summary = {}
for endpoint in cfg["endpoints"]:
    table = bits.filter(pl.col("endpoint") == endpoint)
    n = cfg["endpoint_sizes"][endpoint]
    eligible = table.filter(
        (pl.col("molecules_with_bit") >= cfg["minimum_bit_count"])
        & (pl.col("molecules_with_bit") <= cfg["maximum_frequency"] * n)
    ).with_columns(
        (pl.col("logistic_coefficient_mean").abs() + pl.col("lgbm_oof_mean_abs_contribution"))
        .alias("consensus_magnitude"),
        pl.when(pl.col("logistic_coefficient_mean") > 0)
        .then(pl.lit("positive"))
        .otherwise(pl.lit("negative"))
        .alias("direction"),
    )
    robust = eligible.filter(
        (pl.col("logistic_sign_consistency") == 1.0)
        & (pl.col("enrichment_fold_sign_consistency") >= 0.8)
        & (pl.col("lgbm_gain_fold_count") >= 3)
        & (
            pl.col("logistic_coefficient_mean").sign()
            == pl.col("prevalence_delta").sign()
        )
    )
    for direction in ["positive", "negative"]:
        selected_tables.append(
            robust.filter(pl.col("direction") == direction)
            .sort("consensus_magnitude", descending=True)
            .head(cfg["top_n"])
        )

    logistic_top = set(
        eligible.sort("logistic_coefficient_mean", descending=True)["bit"][:25].to_list()
        + eligible.sort("logistic_coefficient_mean")["bit"][:25].to_list()
    )
    lgbm_top = set(eligible.sort("lgbm_gain_mean", descending=True)["bit"][:50].to_list())
    fdr = eligible.filter(pl.col("fisher_q") <= 0.05).sort("fisher_q")
    endpoint_ablation = ablations.filter(pl.col("endpoint") == endpoint).sort(
        "AUROC", descending=True
    )
    summary[endpoint] = {
        "eligible_bits": eligible.height,
        "robust_consensus_bits": robust.height,
        "fdr_0_05_nonubiquitous_bits": fdr.height,
        "fdr_bits": fdr["bit"].to_list(),
        "logistic_lgbm_top50_overlap": len(logistic_top & lgbm_top),
        "best_ablation_by_auroc": endpoint_ablation.row(0, named=True),
        "all_bits": endpoint_ablation.filter(pl.col("feature_set") == "all_bits").row(
            0, named=True
        ),
    }

selected = pl.concat(selected_tables, how="vertical_relaxed").select(
    "endpoint",
    "direction",
    "bit",
    "molecules_with_bit",
    "positive_rate_with_bit",
    "positive_prevalence",
    "prevalence_delta",
    "odds_ratio",
    "fisher_q",
    "enrichment_fold_sign_consistency",
    "logistic_coefficient_mean",
    "logistic_coefficient_sd",
    "lgbm_gain_mean",
    "lgbm_gain_fold_count",
    "lgbm_oof_mean_contribution_when_present",
    "radius0_occurrences",
    "radius1_occurrences",
    "radius2_occurrences",
    "distinct_fragment_count",
    "top_fragments",
)
selected.write_csv(out / "robust-consensus-bits.csv")
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(selected.write_csv(), flush=True)
print(json.dumps(summary, indent=2), flush=True)
