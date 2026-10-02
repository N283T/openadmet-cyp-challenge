"""Recorded incumbent comparison; promote compact outputs without changing source runs."""

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import polars as pl
from launch_cyp_boltz_incumbent_suite import recorded


def analysis():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    regression = pl.read_parquet(inp / "regression-metrics.parquet")
    tdi = pl.read_csv(inp / "tdi-metrics.csv")
    uncertainty = pl.read_parquet(inp / "paired-uncertainty.parquet")
    summary = json.loads((inp / "regression-summary.json").read_text())
    screens = json.loads((inp / "screen-runs.json").read_text())
    count = 0
    for item in screens:
        fold = pl.read_parquet(inp / f"{item['target']}-fold-metrics.parquet")
        if fold.height != 375:
            raise ValueError("Incomplete regression grid")
        count += fold.height
    if count != 1500 or len(summary["checks"]) != 4:
        raise ValueError("Incomplete incumbent experiment")
    lines = [
        "# Frozen CYP/HEM trunk in the submitted ensembles",
        "",
        "This replaces the old UMAP regression screen for incumbent comparisons. "
        "All four regression readouts were refitted on the exact incumbent parent-series "
        "cohort, labels, confidence bounds and five folds. "
        "Existing trunk extraction and pooling are reused.",
        "",
        "Each incumbent is reconstructed exactly from saved atomic predictions and fold weights. "
        "D6 preserves the submitted two-member core and moment-matched robust branch "
        "at their existing equal top-level weights. Additions and replacements operate "
        "on whole members; no arbitrary fractional blend sweep is used.",
        "",
        "Placement selection and trunk candidate selection exclude the held fold. "
        "The 2C9 bank reruns exhaustive equal-subset selection; 3A4 reruns Caruana "
        "with inner step selection; 1A2 adds one equal member; D6 also evaluates Caruana "
        "over its two preserved branches and the selected trunk member. These are conditional "
        "comparisons on saved OOF, not fully nested base-model refits. Other-fold OOF models "
        "may have seen labels from the held fold during their original base training. "
        "Postselection and multiple comparisons limit uncertainty claims.",
        "",
        "## Regression",
        "",
        "Lower ST-RAE is better. Standalone minimum is best-of-grid; use the "
        "outer-selected placement and bank reselection rows to assess ensemble utility.",
        "",
        "| Target | N | Incumbent | Standalone grid minimum | "
        "Outer-selected placement | Bank reselection |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    pooled = regression.filter(pl.col("fold") == -1)
    for target in sorted(pooled["target"].unique()):
        part = pooled.filter(pl.col("target") == target)
        get = lambda name, part=part: part.filter(pl.col("scenario") == name)["ST_RAE"][0]  # noqa: E731
        standalone = (
            part.filter(pl.col("scenario").str.starts_with("standalone::"))
            .sort("ST_RAE")
            .row(0, named=True)
        )
        lines.append(
            f"| {target} | {standalone['n']} | {get('current'):.6f} | "
            f"{standalone['ST_RAE']:.6f} | {get('outer_selected_placement'):.6f} | "
            f"{get('outer_selected_trunk_bank_reselection'):.6f} |"
        )
    lines += [
        "",
        "Same reselection algorithm without trunk provides a separate attribution control:",
        "",
        "| Target | Original bank reselected | Bank plus trunk reselected |",
        "|---|---:|---:|",
    ]
    for target in sorted(pooled["target"].unique()):
        part = pooled.filter(pl.col("target") == target)
        baseline = part.filter(pl.col("scenario") == "original_bank_reselection")["ST_RAE"][0]
        extended = part.filter(
            pl.col("scenario") == "outer_selected_trunk_bank_reselection"
        )["ST_RAE"][0]
        lines.append(f"| {target} | {baseline:.6f} | {extended:.6f} |")
    lines += [
        "",
        "Conditional parent-group bootstrap intervals (2,000 draws, fixed OOF denominator):",
        "",
        "| Target | Scenario | Delta ST-RAE | Lower 95% | Upper 95% |",
        "|---|---|---:|---:|---:|",
    ]
    for row in uncertainty.sort("target", "scenario").iter_rows(named=True):
        lines.append(
            f"| {row['target']} | {row['scenario']} | {row['delta_ST_RAE']:.6f} | "
            f"{row['lower_95']:.6f} | {row['upper_95']:.6f} |"
        )
    lines += [
        "",
        "## TDI diagnostic",
        "",
        "The original 500 fitted models already match the incumbent UMAP folds and labelled "
        "cohorts; no duplicate TDI fitting was performed. Reconstructed D6 core/rank/tail "
        "stages and A4 seven-member rank ensemble are retained. Trunk readouts enter the D6 "
        "selected core or replace its PCA/raw block; A4 accepts a new equal member or replaces "
        "any named incumbent atom. AP-selected placement excludes the held fold. Historical "
        "labelled-fold MCC points are diagnostic only: the full 6,146-molecule trunk pool is "
        "unavailable, so no call-rate or deployment decision follows.",
        "",
        "| Target | Scenario | AP | AUROC |",
        "|---|---|---:|---:|",
    ]
    for target in ["d6", "a4"]:
        for name in ["current", "outer_selected_AP_placement"]:
            row = tdi.filter((pl.col("target") == target) & (pl.col("scenario") == name)).row(
                0, named=True
            )
            lines.append(f"| {target} | {name} | {row['AP']:.6f} | {row['AUROC']:.6f} |")
    lines += [
        "",
        "## Provenance",
        "",
        f"Regression fits completed: {count}. "
        f"Incumbent reconstruction checks: `{json.dumps(summary['checks'])}`.",
        "",
        "All fixed comparisons, keyed predictions and selections remain in their recorded runs. "
        "The previous UMAP scientific results are retained. "
        "No candidate, submission or upload is created.",
    ]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    for name in [
        "regression-metrics.parquet",
        "paired-uncertainty.parquet",
        "tdi-metrics.csv",
        "regression-selections.parquet",
        "d6-trunk-selections.parquet",
        "a4-trunk-selections.parquet",
    ]:
        shutil.copy2(inp / name, out / name)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--recorded":
        analysis()
        return
    root = Path(__file__).resolve().parents[1]
    directory = Path(sys.argv[1]).resolve()
    runs = root / "artifacts/runs"
    regression_id = (directory / "regression-placement-run.txt").read_text().strip()
    tdi_id = (directory / "tdi-placement-run.txt").read_text().strip()
    reg, tdi = (
        runs / "regression" / regression_id / "outputs",
        runs / "classification" / tdi_id / "outputs",
    )
    inputs = {
        "regression-metrics.parquet": reg / "metrics.parquet",
        "regression-summary.json": reg / "summary.json",
        "regression-selections.parquet": reg / "selections.parquet",
        "paired-uncertainty.parquet": reg / "paired-uncertainty.parquet",
        "tdi-metrics.csv": tdi / "metrics.csv",
        "d6-trunk-selections.parquet": tdi / "d6-trunk-selections.parquet",
        "a4-trunk-selections.parquet": tdi / "a4-trunk-selections.parquet",
        "screen-runs.json": directory / "screen-runs.json",
        "plan.json": directory / "plan.json",
    }
    for item in json.loads((directory / "screen-runs.json").read_text()):
        inputs[f"{item['target']}-fold-metrics.parquet"] = (
            runs / "regression" / item["run"] / "outputs/fold-metrics.parquet"
        )
    # The script name accepts --recorded through a tiny launcher configuration.
    if len(sys.argv) > 2:
        identifier = sys.argv[2]
        result = json.loads((runs / "shared" / identifier / "result.json").read_text())
        if result["status"] != "completed":
            raise ValueError("Only a completed recorded report can be finalized")
    else:
        identifier = recorded(
            directory,
            f"incumbent-comparison-report-{regression_id}",
            "shared",
            {"regression_run": regression_id, "tdi_run": tdi_id},
            inputs,
            "report_cyp_boltz_incumbent_recorded.py",
        )
    (directory / "report-run.txt").write_text(identifier + "\n")
    source = runs / "shared" / identifier / "outputs"
    tracked = (
        root / "artifacts/tracked"
        / f"cyp-boltz-incumbent-{directory.name}-{regression_id.rsplit('-', 1)[-1]}"
    )
    tracked.mkdir(exist_ok=True)
    inventory = []
    for path in source.iterdir():
        if path.is_file() and path.stat().st_size < 50_000_000:
            destination = tracked / path.name
            if destination.exists():
                if destination.read_bytes() != path.read_bytes():
                    raise ValueError("Retained tracked output differs from recorded report")
            else:
                shutil.copy2(path, destination)
            inventory.append(
                {
                    "file": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
            )
    (tracked / "provenance.json").write_text(
        json.dumps(
            {
                "report_run": identifier,
                "regression_run": regression_id,
                "tdi_run": tdi_id,
                "suite": str(directory.relative_to(root)),
                "files": inventory,
            },
            indent=2,
        )
        + "\n"
    )
    shutil.copy2(source / "report.md", directory / "report.md")
    with (root / "docs/research/2026-10.md").open("a") as stream:
        stream.write(
            f"\n### Incumbent-matched trunk redo completed\n\nSuite `{directory.name}` "
            "refits all 1,500 regression fits on incumbent parent cohorts and evaluates "
            "whole-member ensemble placement. TDI reuses matching completed fits. "
            f"Recorded report: `shared/{identifier}`; regression placement: `{regression_id}`; "
            f"TDI placement: `{tdi_id}`. "
            f"[Comparison report](../../{tracked.relative_to(root)}/report.md). "
            "No candidate or upload.\n"
        )


if __name__ == "__main__":
    main()
