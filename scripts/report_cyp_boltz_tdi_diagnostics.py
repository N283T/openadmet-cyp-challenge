"""Render the retained partial-pool diagnostic tables without recomputing experiments."""

import json
import os
from pathlib import Path

import polars as pl

VIEW_DESCRIPTIONS = {
    "ligand_only": "Ligand s mean; ligand-internal z mean/max",
    "pxr_protein": "Protein s mean + ligand s mean; whole-protein-to-ligand z mean/max",
    "pxr_core": "Protein s mean + ligand s mean; core-pocket-to-ligand z mean/max",
    "pxr_extended": "Protein s mean + ligand s mean; extended-pocket-to-ligand z mean/max",
    "core_mixed_hem": "Protein/ligand s means; core pocket and HEM pooled together for ligand z",
    "extended_mixed_hem": "Protein/ligand s means; extended pocket and HEM pooled together",
    "core_separate_hem": (
        "Protein/ligand s means; core-pocket-to-ligand and HEM-to-ligand z separately"
    ),
    "extended_separate_hem": "Extended-pocket-to-ligand and HEM-to-ligand z separately",
    "core_separate_fe": "Core pocket, Fe and non-Fe HEM-to-ligand z separately",
    "core_bidirectional_hem": "Core pocket and HEM separate; both z directions to/from ligand",
    "core_separate_hem_s": "Core pocket/HEM separate z blocks, plus HEM s mean",
    "core_distribution_hem": "Core pocket/HEM separate z mean/std/q10/q90; ligand s mean/std",
}


def table(frame, columns):
    def cell(value):
        if value is None:
            return "—"
        if isinstance(value, float):
            return f"{value:.6f}"
        return str(value).replace("|", "/").replace("\n", " ")

    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    lines.extend(
        "| " + " | ".join(cell(row.get(col)) for col in columns) + " |" for row in frame.to_dicts()
    )
    return "\n".join(lines)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    summary = json.loads((inp / "summary.json").read_text())
    metrics = pl.read_parquet(inp / "metrics.parquet")
    selected = metrics.filter(
        (pl.col("scenario") == "current")
        | (
            pl.col("scenario").str.contains("selected|single_joint_rate")
            & ~pl.col("scenario").str.contains("::fold")
        )
    )
    intervals = pl.read_parquet(inp / "intervals.parquet")
    choices = pl.read_parquet(inp / "selections.parquet")
    correlations = pl.read_parquet(inp / "incumbent-correlations.parquet")
    pairs = pl.read_parquet(inp / "trunk-correlations.parquet")
    standalone = metrics.filter(pl.col("scenario").str.starts_with("standalone::"))
    leading = pl.concat(
        [
            standalone.filter(pl.col("target") == t).sort("AP", descending=True).head(5)
            for t in ["d6", "a4"]
        ]
    )
    lines = [
        "# Frozen CYP/HEM trunk: labelled-only TDI integration diagnostics",
        "",
        f"Diagnostic source run: `{cfg['source_run']}`. No new model fits or trunk extractions.",
        "",
        "The user explicitly deferred full-pool extraction because of cost, effort and storage. "
        "This analysis looks for signals worth completing the operating pool later. It cannot "
        "set a deployment rate or justify a submitted call vector. Regression is unchanged.",
        "",
        "## Cohort and controls",
        "",
        table(
            pl.DataFrame(summary["coverage"]),
            [
                "target",
                "labelled_rows",
                "total_pool_rows",
                "missing_rows",
                "reconstruction_max_diff",
            ],
        ),
        "",
        "Each trunk recipe covers exactly the incumbent labelled molecule IDs and UMAP folds. "
        "Measurement labels and saved fold assignments are checked explicitly. The exact saved "
        "incumbent composition is reproduced. Missing classification labels are never negatives.",
        "",
        "D6 preserves its saved core-member choices, PCA/raw TabPFN block, rank transform and "
        "Morgan logistic/SVM/complement-NB tail stages. A4 preserves its seven rank members: "
        "historical stacked TabPFN with Mordred, pretrained molecular 600-dimensional embeddings "
        "with PCA256/raw600 and TabICL/TabPFN, external-CYP classification embeddings plus Mordred "
        "with LightGBM, and Morgan plus PXR/ADMET auxiliary predictions with LightGBM.",
        "",
        "The no-trunk control uses the same saved original scores and tests adding an existing "
        "original atom as a core/equal-mean member; it diagnoses weight/composition changes. "
        "It is not an unrestricted reselection of every historical model bank. Molecular-only "
        "CheMeleon2048 + Mordred1613 readouts are separately included as ancillary placements.",
        "",
        "## Conditional selection procedures",
        "",
        "Single addition/replacement is chosen on the other four saved OOF folds, separately by "
        "AP and target-prior MCC at the historical diagnostic rate. The joint policy selects "
        "both single placement and rate on those four folds. All include keeping the incumbent.",
        "",
        "Multi-block selection uses inner held folds among those four to select block sizes "
        "1/2/3/5 and either one averaged member, multiple individual additions, or 10/25/50% "
        "mixing of incumbent ranks with trunk-block ranks. Distinct readouts are greedily "
        "selected on inner training folds, then reselected on all four folds for the outer score. "
        "Greedy replacement preserves existing member positions and stops when training-fold "
        "quality no longer improves: D6 replaces only its PCA/raw block; A4 can replace up to "
        "three distinct members with distinct readouts.",
        "",
        "All methods are conditional saved-OOF analyses, not fully nested base-model refits. "
        "Other-fold OOF base models may have used evaluated-fold labels "
        "in their original training. "
        "Multiple trunks here are readouts/views of the same frozen extraction, not independently "
        "computed Boltz trunks or receptor conformations.",
        "",
        "## Selected policies",
        "",
        table(
            selected,
            [
                "target",
                "scenario",
                "AP",
                "AUROC",
                "MCC_target_prior",
                "added",
                "removed",
                "positive_added",
                "positive_removed",
            ],
        ),
        "",
        "All calls above rank labelled molecules only. Added/removed counts are diagnostic "
        "training calls, not blind-test changes. Raw probability Brier scores, fold-wise metrics, "
        "all fixed single placements and prior/rate curves are in the accompanying tables.",
        "",
        "## Standalone representation screen",
        "",
        table(leading, ["target", "scenario", "AP", "AUROC", "MCC_target_prior"]),
        "",
        "These are global best-of-grid standalone readouts, shown only for representation "
        "diagnosis. They are not estimates of selected ensemble improvement.",
        "",
        "## Paired conditional uncertainty",
        "",
        table(
            intervals,
            [
                "target",
                "scenario",
                "reference",
                "delta_AP",
                "AP_lower",
                "AP_upper",
                "delta_MCC",
                "MCC_lower",
                "MCC_upper",
            ],
        ),
        "",
        "2,000 paired UMAP-group resamples retain selected scores and calls. They do not rerank, "
        "refit selectors or account for searching the grid. Intervals are conditional descriptive "
        "uncertainty, not independent confirmation or multiplicity-adjusted evidence.",
        "",
        "## Fold choices",
        "",
        table(choices, ["target", "objective", "family", "fold", "recipe", "rate"]),
        "",
        "## Prediction and error complementarity",
        "",
    ]
    for target in ["d6", "a4"]:
        c = correlations.filter(
            (pl.col("target") == target)
            & (pl.col("reference") == "incumbent")
            & (pl.col("fold") == -1)
        )
        top = c.sort("prediction_spearman").head(8)
        lines.extend(
            [
                f"### {target}: least-correlated readouts versus incumbent",
                "",
                table(
                    top,
                    [
                        "model",
                        "prediction_spearman",
                        "rank_residual_pearson",
                        "rank_squared_loss_pearson",
                    ],
                ),
                "",
            ]
        )
        p = pairs.filter(pl.col("target") == target)
        lines.append(
            f"Trunk pair prediction Spearman range: {p['prediction_spearman'].min():.6f} "
            f"to {p['prediction_spearman'].max():.6f}; "
            f"median {p['prediction_spearman'].median():.6f}."
        )
        lines.append("")
    lines.extend(
        [
            "Residuals here are within-fold percentile ranks minus binary "
            "labels. Their shared label "
            "term can produce high correlation, so low prediction correlation "
            "alone is not evidence "
            "of useful error diversity. Squared-rank-loss correlations and held-fold call swaps "
            "provide complementary diagnostics. These are not calibrated-probability residuals.",
            "",
            "## Feature interpretation",
            "",
            table(
                pl.DataFrame(
                    [
                        {"view": view, "features": description}
                        for view, description in VIEW_DESCRIPTIONS.items()
                    ]
                ),
                ["view", "features"],
            ),
            "",
            "`trunk_molecular` additionally uses frozen CheMeleon2048 + Mordred1613. "
            "`lgbm_binary` is the fixed LightGBM classifier; `tabpfn_top500` is TabPFN 2.6 "
            "with eight estimators and training-fold-only LightGBM gain "
            "selection up to 500 columns. "
            "Boltz z statistics are latent embeddings, not measured distances.",
            "",
            "## Decision limits",
            "",
            "A useful follow-up signal should improve the actual incumbent, survive held-fold "
            "selection and no-trunk composition controls, retain useful "
            "labelled positive additions "
            "at nearby rates, and show more than a global grid maximum. Fold "
            "instability, rate-grid "
            "boundary selection and wide paired intervals lower the case for storage-expensive "
            "full-pool extraction. AP increases alone do not establish MCC or "
            "leaderboard improvement.",
            "",
            "The historical D6 AP-versus-control replay and A4 249-to-275 rate probe both failed "
            "to reproduce prospective leaderboard directions. Thus even an encouraging signal "
            "here is only grounds to consider completing the pool, not submission authorization.",
            "",
        ]
    )
    (out / "report.md").write_text("\n".join(lines))
    (out / "summary.json").write_text(
        json.dumps(
            {"source_run": cfg["source_run"], "deployment_decision_supported": False}, indent=2
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
