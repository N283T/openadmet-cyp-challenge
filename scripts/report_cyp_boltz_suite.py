"""Record paired OOF analysis of the full fixed CYP/HEM pooling comparison."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]


def analyse():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / "measurements.parquet")
    text = [
        "# CYP frozen Boltz trunk pooling screen",
        "",
        "This is an exploratory fixed-grid comparison, not a submission decision. "
        "The Boltz trunk is frozen; all tensors include HEM. 'No HEM' variants omit "
        "HEM from readout pooling, not from the trunk input.",
        "",
        "All models use the original five UMAP folds. Molecular controls use "
        "CheMeleon2048 + Mordred1613. Imputation and TabPFN top500 gain selection "
        "are fitted only on each training fold. Test predictions average five models.",
        "",
        "Trunk features contain no predicted coordinates. Pair features are latent "
        "representations and are not interpreted as measured contacts.",
        "",
        "Rankings and blend weights below are exploratory. No nested recipe selection "
        "or independent confirmation is claimed. Paired intervals are conditional "
        "on these fitted models, without multiplicity correction.",
        "",
    ]
    all_metrics, uncertainty = [], []
    for index, spec in enumerate(cfg["screens"]):
        prefix = f"screen-{index}"
        metrics_path = inp / f"{prefix}-aggregate.parquet"
        text += [
            f"## {spec['target']} / {spec['track']}",
            "",
            f"Run: `{spec['run']}`; launcher exit code: {spec['exit_code']}.",
            "",
        ]
        if spec.get("recovery_run"):
            text += [f"Export recovered without refitting in `{spec['recovery_run']}`.", ""]
        if not metrics_path.exists():
            text += [
                "No complete recipe metrics were available. Inspect retained run failure outputs.",
                "",
            ]
            continue
        metrics = pl.read_parquet(metrics_path)
        scores = pl.read_parquet(inp / f"{prefix}-folds.parquet")
        summary = (
            json.loads((inp / f"{prefix}-summary.json").read_text())
            if (inp / f"{prefix}-summary.json").exists()
            else {}
        )
        text += [
            (
                f"Completed recipes: {metrics.height}. "
                f"Recorded failures: {len(summary.get('failures', []))}."
            ),
            "",
        ]
        if spec["track"] == "classification":
            text += [
                "TDI is a labelled-only diagnostic within the target-extracted pool. "
                "The complete 6,146-molecule operating pool is absent, so AP/AUROC "
                "and MCC at 0.5 do not support call-rate selection or deployment.",
                "",
            ]
            ranking = metrics.sort("AP", descending=True)
            text += ["| Recipe | AP | AUROC | MCC at 0.5 |", "|---|---:|---:|---:|"]
            for r in ranking.iter_rows(named=True):
                text.append(
                    f"| {r['model']} | {r['AP']:.4f} | {r['AUROC']:.4f} | {r['MCC_at_0.5']:.4f} |"
                )
        else:
            ranking = metrics.sort("ST_RAE")
            text += [
                "| Recipe | ST-RAE | Spearman | Blend 10% | Blend 25% | Blend 50% |",
                "|---|---:|---:|---:|---:|---:|",
            ]
            for r in ranking.iter_rows(named=True):
                blends = " | ".join(
                    f"{r[f'blend_{w}_ST_RAE']:.4f}"
                    if r.get(f"blend_{w}_ST_RAE") is not None
                    else "—"
                    for w in [0.1, 0.25, 0.5]
                )
                text.append(
                    f"| {r['model']} | {r['ST_RAE']:.4f} | {r['Spearman']:.4f} | {blends} |"
                )
            if not (inp / f"{prefix}-predictions.parquet").exists():
                text += [
                    "",
                    "The run stopped before its combined prediction bank was saved; "
                    "partial metrics are retained without paired uncertainty.",
                    "",
                ]
                all_metrics.append(
                    metrics.with_columns(
                        pl.lit(spec["target"]).alias("target"), pl.lit(spec["track"]).alias("track")
                    )
                )
                scores.write_parquet(out / f"{prefix}-fold-metrics.parquet")
                continue
            bank = (
                pl.read_parquet(inp / f"{prefix}-predictions.parquet")
                .filter(pl.col("role") == "train")
                .join(
                    labels.select("compound_id", "endpoint", "value", "conf_low", "conf_high"),
                    on=["compound_id", "endpoint"],
                    how="inner",
                    validate="m:1",
                )
            )
            text += [
                "",
                "Paired group-bootstrap differences versus the same-readout molecular control "
                "are saved in `paired-uncertainty.parquet`. Negative ST-RAE differences favour "
                "the trunk recipe. The reference denominator is fixed at the complete OOF "
                "label mean; intervals describe grouped error differences at that denominator.",
                "",
            ]
            for r in metrics.iter_rows(named=True):
                if r["mode"] == "molecular":
                    continue
                control_name = f"molecular_only--molecular--{r['method']}"
                control = bank.filter(pl.col("model") == control_name).select(
                    "compound_id", "fold", pl.col("prediction").alias("control")
                )
                paired = bank.filter(pl.col("model") == r["model"]).join(
                    control, on=["compound_id", "fold"], validate="1:1"
                )
                y, low, high = (paired[c].to_numpy() for c in ["value", "conf_low", "conf_high"])
                p, c = paired["prediction"].to_numpy(), paired["control"].to_numpy()

                def err(prediction, low=low, high=high):
                    return np.maximum(low - prediction, 0) + np.maximum(prediction - high, 0)

                denominator = err(y.mean()).sum()
                for weight in [1.0, 0.1, 0.25, 0.5]:
                    delta = err((1 - weight) * c + weight * p) - err(c)
                    groups = (
                        paired.select("group")
                        .with_columns(pl.Series("delta", delta))
                        .group_by("group")
                        .agg(pl.col("delta").sum())["delta"]
                        .to_numpy()
                    )
                    rng = np.random.default_rng(20261001)
                    boot = np.array(
                        [
                            rng.choice(groups, len(groups), replace=True).sum() / denominator
                            for _ in range(2000)
                        ]
                    )
                    by_fold = (
                        paired.select("fold")
                        .with_columns(pl.Series("delta", delta))
                        .group_by("fold")
                        .agg(pl.col("delta").sum())
                    )
                    uncertainty.append(
                        {
                            "target": spec["target"],
                            "model": r["model"],
                            "weight": weight,
                            "delta_ST_RAE": float(delta.sum() / denominator),
                            "lower_95": float(np.quantile(boot, 0.025)),
                            "upper_95": float(np.quantile(boot, 0.975)),
                            "improving_folds": int((by_fold["delta"] < 0).sum()),
                            "groups": len(groups),
                            "conditional_exploratory": True,
                        }
                    )
        all_metrics.append(
            metrics.with_columns(
                pl.lit(spec["target"]).alias("target"), pl.lit(spec["track"]).alias("track")
            )
        )
        text += [
            "",
            "All fixed recipes are shown; the minimum/maximum is not an unbiased estimate "
            "of a selected model's future performance.",
            "",
        ]
        scores.with_columns(pl.lit(spec["target"]).alias("target")).write_parquet(
            out / f"{prefix}-fold-metrics.parquet"
        )
    if all_metrics:
        pl.concat(all_metrics, how="diagonal_relaxed").write_parquet(out / "all-metrics.parquet")
    if uncertainty:
        pl.DataFrame(uncertainty).write_parquet(out / "paired-uncertainty.parquet")
    text += [
        "## Failures and scope",
        "",
        json.dumps(cfg["failures"], indent=2),
        "",
        "Original large tensor and pooled feature payloads stay local with recorded input hashes. "
        "Completed runs and negative results are retained. No candidate or submission is uploaded.",
        "",
    ]
    (out / "report.md").write_text("\n".join(text))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", nargs="?", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.execute:
        analyse()
        return
    directory = args.directory.resolve()
    screens = (
        json.loads((directory / "screen-runs.json").read_text())
        if (directory / "screen-runs.json").exists()
        else []
    )
    # The fitting source snapshot cannot be changed during a running screen.
    # Recover known export-only errors from saved shards before reporting.
    for spec in screens:
        if spec["exit_code"] == 0 or not spec["run"]:
            continue
        sidecar = directory / "recovery-runs.json"
        recovered = json.loads(sidecar.read_text()) if sidecar.exists() else {}
        source = ROOT / "artifacts/runs" / spec["track"] / spec["run"]
        console = source / "console.log"
        if (
            spec["run"] not in recovered
            and console.exists()
            and ("type Float32 is incompatible with expected type Float64" in console.read_text())
        ):
            subprocess.run(
                [
                    "pixi",
                    "run",
                    "python",
                    "scripts/recover_cyp_boltz_predictions.py",
                    spec["track"],
                    spec["run"],
                    "--suite",
                    str(directory),
                ],
                cwd=ROOT,
                check=True,
            )
    sidecar = directory / "recovery-runs.json"
    recovered = json.loads(sidecar.read_text()) if sidecar.exists() else {}
    source_failures = json.loads((directory / "failures.json").read_text())
    resolved_pairs = {(s["target"], s["track"]) for s in screens if s["run"] in recovered}
    failures = [
        f for f in source_failures if (f.get("target"), f.get("track")) not in resolved_pairs
    ]
    files = {
        "suite-plan.json": directory / "plan.json",
        "suite-failures.json": directory / "failures.json",
        "measurements.parquet": ROOT
        / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    }
    if sidecar.exists():
        files["suite-recoveries.json"] = sidecar
    for i, spec in enumerate(screens):
        if not spec["run"]:
            continue
        source_id = spec["run"]
        if source_id in recovered:
            recovery = recovered[source_id]
            result = json.loads(
                (
                    ROOT / "artifacts/runs" / recovery["track"] / recovery["run"] / "result.json"
                ).read_text()
            )
            if result["status"] != "completed":
                raise ValueError("Recovery is not complete")
            spec["recovery_run"] = recovery["run"]
            files[f"screen-{i}-recovery-result.json"] = (
                ROOT / "artifacts/runs" / recovery["track"] / recovery["run"] / "result.json"
            )
        output = (
            ROOT
            / "artifacts/runs"
            / spec["track"]
            / spec.get("recovery_run", source_id)
            / "outputs"
        )
        for name, alias in [
            ("aggregate-metrics.parquet", "aggregate.parquet"),
            ("fold-metrics.parquet", "folds.parquet"),
            ("predictions.parquet", "predictions.parquet"),
            ("summary.json", "summary.json"),
        ]:
            if (output / name).exists():
                files[f"screen-{i}-{alias}"] = output / name
    config = ROOT / "configs/boltz-trunk" / f"report-{directory.name}.json"
    config.write_text(
        json.dumps(
            {"screens": screens, "failures": failures, "source_failures": source_failures}, indent=2
        )
        + "\n"
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Analyse all fixed CYP/HEM pooling ablations with paired grouped OOF comparisons",
        "--config",
        str(config),
    ]
    for alias, path in files.items():
        command += ["--input", f"{alias}={path}"]
    result = subprocess.run(
        command + ["--", "python", "scripts/report_cyp_boltz_suite.py", "--execute"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    (directory / "report-launch.log").write_text(result.stdout + result.stderr)
    result.check_returncode()
    identifier = next(
        line.split("=", 1)[1] for line in result.stdout.splitlines() if line.startswith("run_id=")
    )
    (directory / "report-run.txt").write_text(identifier + "\n")
    report = ROOT / "artifacts/runs/shared" / identifier / "outputs/report.md"
    shutil.copy2(report, directory / "report.md")
    tracked = ROOT / "artifacts/tracked" / f"cyp-boltz-pooling-{directory.name}"
    tracked.mkdir()
    from cyp_challenge.storage import digest

    promoted = []
    for path in sorted(report.parent.iterdir()):
        if path.is_file() and path.stat().st_size < 50_000_000:
            shutil.copy2(path, tracked / path.name)
            promoted.append(
                {"file": path.name, "sha256": digest(path), "bytes": path.stat().st_size}
            )
    (tracked / "provenance.json").write_text(
        json.dumps(
            {
                "report_run": identifier,
                "suite": str(directory.relative_to(ROOT)),
                "screens": screens,
                "files": promoted,
                "large_features_and_model_assets": (
                    "Retained in source runs; not copied into ordinary Git"
                ),
            },
            indent=2,
        )
        + "\n"
    )
    log = ROOT / "docs/research/2026-10.md"
    completed_screens = sum(s["exit_code"] == 0 or bool(s.get("recovery_run")) for s in screens)
    with log.open("a") as stream:
        stream.write(
            f"\n### Completion of suite `{directory.name}`\n\n"
            f"Recorded report run: `shared/{identifier}`. "
            f"[Full fixed-grid comparison](../../artifacts/tracked/{tracked.name}/report.md). "
            f"Completed screens: {completed_screens}/{len(screens)}. "
            "All results, including failures and negative comparisons, are retained. "
            "No candidate or submission was created.\n"
        )
    print(f"report_run={identifier}\nreport={directory / 'report.md'}", flush=True)


if __name__ == "__main__":
    main()
