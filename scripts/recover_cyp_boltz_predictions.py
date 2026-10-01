"""Recover fully fitted screen prediction shards in a new immutable recorded run."""

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]


def combine_predictions(frames):
    """Normalize prediction precision without changing numerical values."""
    return pl.concat([f.with_columns(pl.col("prediction").cast(pl.Float64)) for f in frames])


def execute():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    frames = [
        pl.read_parquet(inp / name).with_columns(pl.col("prediction").cast(pl.Float64))
        for name in cfg["shards"]
    ]
    bank = combine_predictions(frames)
    keys = ["compound_id", "endpoint", "split_id", "fold", "model", "role"]
    if bank.select(keys).unique().height != bank.height:
        raise ValueError("Duplicate prediction keys")
    if not np.isfinite(bank["prediction"].to_numpy()).all():
        raise ValueError("Nonfinite predictions")
    counts = bank.filter(pl.col("role") == "test").group_by("model", "fold").len()
    if counts.height != cfg["expected_fits"] or not (counts["len"] == 750).all():
        raise ValueError("Incomplete model/fold test coverage")
    if set(bank["fold"]) != set(range(5)):
        raise ValueError("Missing folds")
    recipes = pl.read_parquet(inp / "aggregate.parquet")
    fold_metrics = pl.read_parquet(inp / "fold-metrics.parquet")
    if recipes.height != cfg["expected_recipes"] or fold_metrics.height != cfg["expected_fits"]:
        raise ValueError("Not all fits completed; do not claim full export recovery")
    bank.write_parquet(out / "predictions.parquet")
    recipes.write_parquet(out / "aggregate-metrics.parquet")
    fold_metrics.write_parquet(out / "fold-metrics.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "source_run": cfg["source_run"],
                "endpoint": bank["endpoint"][0],
                "successful_recipes": recipes.height,
                "planned_recipes": cfg["expected_recipes"],
                "failures": [],
                "source_export_error": "Float32/Float64 concat mismatch",
                "recovery": (
                    "Cast saved prediction shards to Float64; no refitting or prediction changes"
                ),
                "exploratory_only": True,
                "classification_deployment_supported": False,
                "prediction_rows": bank.height,
                "fits": fold_metrics.height,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        (
            f"Recovered {recipes.height} recipes / {fold_metrics.height} fits / "
            f"{bank.height} predictions"
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("track", nargs="?", choices=["regression", "classification"])
    parser.add_argument("source_run", nargs="?")
    parser.add_argument("--suite", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.execute:
        execute()
        return
    run = ROOT / "artifacts/runs" / args.track / args.source_run
    result = json.loads((run / "result.json").read_text())
    if (
        result["status"] != "failed"
        or "type Float32 is incompatible with expected type Float64"
        not in (run / "console.log").read_text()
    ):
        raise ValueError("Recovery is restricted to the known terminal concat error")
    shards = sorted((run / "outputs/predictions").glob("*--fold*.parquet"))
    fits = 375 if args.track == "regression" else 250
    if len(shards) != fits:
        raise ValueError(f"Expected {fits} completed prediction shards, got {len(shards)}")
    config = ROOT / "configs/boltz-trunk" / f"recover-{args.source_run}.json"
    aliases = {f"shard-{i:04d}.parquet": p for i, p in enumerate(shards)}
    payload = {
        "source_run": args.source_run,
        "track": args.track,
        "shards": list(aliases),
        "expected_fits": fits,
        "expected_recipes": 75 if args.track == "regression" else 50,
    }
    if not config.exists():
        config.write_text(json.dumps(payload, indent=2) + "\n")
    aliases.update(
        {
            "aggregate.parquet": run / "outputs/aggregate-metrics.parquet",
            "fold-metrics.parquet": run / "outputs/fold-metrics.parquet",
            "source-result.json": run / "result.json",
            "source-config.json": run / "config.json",
            "source-manifest.json": run / "manifest.json",
        }
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        args.track,
        "--purpose",
        f"Recover complete frozen trunk prediction bank from {args.source_run} without refitting",
        "--config",
        str(config),
    ]
    for alias, path in aliases.items():
        command += ["--input", f"{alias}={path}"]
    completed = subprocess.run(
        command + ["--", "python", "scripts/recover_cyp_boltz_predictions.py", "--execute"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    print(completed.stdout[-1500:], completed.stderr[-1500:], flush=True)
    completed.check_returncode()
    identifier = next(
        line.split("=", 1)[1]
        for line in completed.stdout.splitlines()
        if line.startswith("run_id=")
    )
    if args.suite:
        sidecar = args.suite / "recovery-runs.json"
        existing = json.loads(sidecar.read_text()) if sidecar.exists() else {}
        existing[args.source_run] = {"track": args.track, "run": identifier}
        sidecar.write_text(json.dumps(existing, indent=2) + "\n")
    print(f"recovery_run={identifier}", flush=True)


if __name__ == "__main__":
    main()
