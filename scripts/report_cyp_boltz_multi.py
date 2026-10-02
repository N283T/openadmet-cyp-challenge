"""Record and retain a compact report of conditional multiple-trunk placement."""

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import polars as pl
from launch_cyp_boltz_incumbent_suite import recorded


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    metrics = pl.read_parquet(inp / "metrics.parquet").filter(pl.col("fold") == -1)
    pairs = pl.read_parquet(inp / "trunk-pair-correlations.parquet").filter(pl.col("fold") == -1)
    intervals = pl.read_parquet(inp / "paired-uncertainty.parquet")
    selections = pl.read_parquet(inp / "selections.parquet")
    prior = pl.read_parquet(inp / "prior-metrics.parquet").filter(pl.col("fold") == -1)
    lines = [
        "# Multiple trunk models in the incumbent ensemble",
        "",
        f"Analysis run: `regression/{cfg['analysis_run']}`. "
        "No base-model refitting, candidate or upload.",
        "",
        "All four submitted regression cohorts, labels, confidence bounds and parent-series folds "
        "are reused exactly. Each target has 72 trunk recipes: 12 views, with or without molecular "
        "features, and three readouts. Pairwise correlations describe pooled OOF.",
        "",
        "## Trunk diversity",
        "",
        "| Target | Prediction correlation median (10th–90th percentile) "
        "| Residual correlation median |",
        "|---|---:|---:|",
    ]
    for target in cfg["targets"]:
        f = pairs.filter(pl.col("target") == target)
        a, b = f["prediction_pearson"], f["residual_pearson"]
        lines.append(
            f"| {target} | {a.median():.3f} ({a.quantile(0.1):.3f}–{a.quantile(0.9):.3f}) "
            f"| {b.median():.3f} |"
        )
    lines.extend(
        [
            "",
            "## Held-fold results",
            "",
            "ST_RAE is lower-is-better. Counts and blending strength use inner held-fold "
            "error within the other four folds. Greedy block fitting uses distinct trunk recipes, "
            "with block sizes 1, 2, 3, 5 and 10. A block is one prediction column, blended "
            "at weight 0, .05, .1, .2, .3, .5 or 1. Multiple additions give each recipe one "
            "equal-strength slot and preserve incumbent relative weights. Replacement greedily "
            "replaces up to three distinct incumbent atoms at their saved weights; its cap is "
            "inner-fold selected. D6 replaces at most one branch to retain its hierarchy. "
            "Full-bank selection uses Caruana over incumbent atoms plus all 72 trunks, with steps "
            "5, 10, 20 or 50 chosen by inner folds; an incumbent-only control is included.",
            "",
            "| Target | Current | Previous one-member placement "
            "| Previous one-trunk bank extension | Trunk block added | Multiple added "
            "| Multiple replaced | Full bank | Full-bank control |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for target in cfg["targets"]:
        m = {
            r["scenario"]: r["ST_RAE"]
            for r in metrics.filter(pl.col("target") == target).to_dicts()
        }
        p = {
            r["scenario"]: r["ST_RAE"] for r in prior.filter(pl.col("target") == target).to_dicts()
        }
        values = [
            m["current"],
            p["outer_selected_placement"],
            p["outer_selected_trunk_bank_reselection"],
            m["trunk_block_added"],
            m["multiple_added"],
            m["multiple_replaced"],
            m["full_bank"],
            m["original_bank_control"],
        ]
        lines.append(f"| {target} | " + " | ".join(f"{v:.6f}" for v in values) + " |")
    lines.extend(
        [
            "",
            "## Paired uncertainty",
            "",
            "Parent-series group bootstrap, 2,000 draws; conditional on saved OOF and "
            "fitted selections. Intervals do not cover repeating selection or choosing "
            "the best reported method. All improving methods have intervals spanning zero.",
            "",
            "| Target | Method | Reference | Delta ST_RAE | Conditional 95% interval |",
            "|---|---|---|---:|---:|",
        ]
    )
    for row in intervals.filter(pl.col("scenario") != "trunk_block").to_dicts():
        lines.append(
            f"| {row['target']} | {row['scenario']} | {row['reference']} "
            f"| {row['delta_ST_RAE']:.6f} "
            f"| [{row['lower_95']:.6f}, {row['upper_95']:.6f}] |"
        )
    lines.extend(
        [
            "",
            "## Selected compositions",
            "",
            "Recipe identifiers below are exact reproducibility keys. In the 1A2 replacements, "
            "the main pair is ligand-only TabPFN and pocket/HEM bidirectional TabPFN; "
            "some folds also add pocket/HEM bidirectional LGBM. No molecular features are "
            "used in these replacement recipes.",
            "",
        ]
    )
    for row in selections.to_dicts():
        lines.append(
            f"- {row['target']}, fold {row['fold']}, {row['method']}, "
            f"parameter {row['parameter']}: `{row['choices']}`"
        )
    lines.extend(
        [
            "",
            "## Interpretation and limits",
            "",
            "1A2 multiple replacement gives 0.665644 versus current 0.670833 and previous "
            "single replacement 0.669032. Trunk aggregation gives smaller gains. D6 full-bank "
            "selection gives 0.815874, comparable to the previous single-trunk extension "
            "0.815745; extra trunks do not establish further improvement. C9 changes are "
            "negligible or adverse; A4 has no benefit. Trunk-only aggregation remains worse "
            "than the incumbent. Correlated errors limit the benefit of adding similar views.",
            "",
            "This is conditional saved-OOF selection, not fully nested base-model refitting. "
            "Base predictions in other folds can have been trained using held-fold labels. "
            "No test-set prediction, deployment choice or submission was produced.",
        ]
    )
    (out / "report.md").write_text("\n".join(lines) + "\n")


def launch():
    root = Path(__file__).resolve().parents[1]
    suite = root / "artifacts/boltz-trunk/experiments/20261001T065518Z-incumbent"
    identifier = (suite / "multi-trunk-run.txt").read_text().strip()
    source = root / "artifacts/runs/regression" / identifier
    prior = (suite / "regression-placement-run.txt").read_text().strip()
    cfg = {
        "analysis_run": identifier,
        "prior_run": prior,
        "targets": ["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"],
    }
    inputs = {
        name: source / "outputs" / name
        for name in [
            "metrics.parquet",
            "selections.parquet",
            "paired-uncertainty.parquet",
            "trunk-pair-correlations.parquet",
        ]
    }
    inputs["prior-metrics.parquet"] = (
        root / "artifacts/runs/regression" / prior / "outputs/metrics.parquet"
    )
    run = recorded(suite, "multi-trunk-report", "regression", cfg, inputs, Path(__file__).name)
    tracked = root / "artifacts/tracked" / f"cyp-boltz-multi-report-{run}"
    tracked.mkdir()
    report = root / "artifacts/runs/regression" / run / "outputs/report.md"
    shutil.copy2(report, tracked / "report.md")
    (tracked / "provenance.json").write_text(
        json.dumps(
            {
                "run": run,
                "analysis_run": identifier,
                "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"report={tracked / 'report.md'}")


if __name__ == "__main__":
    launch() if "--launch" in sys.argv else main()
