"""Redo regression on submitted parent cohorts and evaluate whole-member placement."""

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs"
OLD = ROOT / "artifacts/boltz-trunk/experiments/20260930T232717Z"
BANKS = {
    "cyp1a2": (
        "20260918T034159381060Z-46314b50",
        "outputs/1A2/add_dense_and_fast_pool",
        "equal_all",
        1372,
    ),
    "cyp2c9": ("20260917T123601175420Z-e1560701", "outputs", "equal_subset", 1255),
    "cyp3a4": ("20260917T130718808881Z-dfa02730", "outputs", "caruana", 2300),
}


def recorded(directory, name, track, config, inputs, script):
    config_path = ROOT / "configs/boltz-trunk" / f"{name}-{directory.name}.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        track,
        "--purpose",
        name,
        "--config",
        str(config_path),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    log = directory / f"{name}.log"
    with log.open("x") as stream:
        code = subprocess.run(
            command + ["--", "python", f"scripts/{script}"],
            cwd=ROOT,
            stdout=stream,
            stderr=subprocess.STDOUT,
        ).returncode
    ids = [
        line.split("=", 1)[1] for line in log.read_text().splitlines() if line.startswith("run_id=")
    ]
    if code or not ids:
        raise RuntimeError(f"{name} failed; inspect {log}")
    return ids[-1]


def main():
    directory = (
        ROOT
        / "artifacts/boltz-trunk/experiments"
        / (datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ") + "-incumbent")
    )
    directory.mkdir()
    print(f"suite_dir={directory}", flush=True)
    pools = json.loads((OLD / "pool-runs.json").read_text())
    plan = {
        "source_suite": str(OLD.relative_to(ROOT)),
        "regression_fits": 1500,
        "regression_split": "20260911T104703655337Z-a3923255:series",
        "cohort": "Exact submitted OOF labels, bounds, IDs and folds",
        "pool_runs": pools,
        "TDI": (
            "Reuse 500 fits on matching UMAP folds; whole-member placement; "
            "labelled diagnostic only"
        ),
        "selection": "Conditional OOF selection, no fully nested base-model refits",
        "upload": False,
    }
    (directory / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    assignments = pl.read_parquet(
        RUNS / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    )
    for target in ["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"]:
        if target == "cyp2d6":
            bank = pl.read_parquet(
                RUNS / "regression/20260926T085959697226Z-0dbf57d4/outputs/oof.parquet"
            ).filter(
                (pl.col("split_id") == "parent")
                & (pl.col("robust_delta") == 0.5)
                & (pl.col("constraint_weight") == 0.35)
                & (pl.col("placement") == "fold_moment_match")
                & (pl.col("blend_weight") == 0.5)
            )
            count = 1457
        else:
            run, output, method, count = BANKS[target]
            bank = pl.read_parquet(RUNS / "regression" / run / output / "oof.parquet").filter(
                pl.col("method") == method
            )
        features = pl.read_parquet(
            RUNS / "shared" / pools[target] / "outputs/features.parquet",
            columns=["compound_id", "role"],
        )
        joined = bank.join(features, on="compound_id", validate="1:1").join(
            assignments.select("compound_id", pl.col("fold").alias("assignment_fold")),
            on="compound_id",
            validate="1:1",
        )
        if (
            joined.height != count
            or joined.filter(pl.col("fold") != pl.col("assignment_fold")).height
            or joined.filter(pl.col("role") != "train").height
        ):
            raise ValueError(f"Preflight cohort mismatch: {target}")
        print(f"preflight_passed={target}:{count}", flush=True)
    screens = []
    for target in ["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"]:
        log = directory / f"screen-{target}-regression.log"
        with log.open("x") as stream:
            code = subprocess.run(
                [
                    "pixi",
                    "run",
                    "python",
                    "scripts/launch_cyp_boltz_screen.py",
                    target,
                    pools[target],
                    "--protocol",
                    "parent-ensemble",
                ],
                cwd=ROOT,
                stdout=stream,
                stderr=subprocess.STDOUT,
            ).returncode
        ids = [
            line.split("=", 1)[1]
            for line in log.read_text().splitlines()
            if line.startswith("run_id=")
        ]
        screens.append(
            {
                "target": target,
                "track": "regression",
                "run": ids[-1] if ids else None,
                "exit_code": code,
            }
        )
        (directory / "screen-runs.json").write_text(json.dumps(screens, indent=2) + "\n")
        if code:
            raise RuntimeError(f"Screen failed: {target}; {log}")
    inputs = {}
    targets = {}
    for target, (run, output, method, count) in BANKS.items():
        source = RUNS / "regression" / run
        summary = source / output / "summary.json"
        for alias, path in [
            (f"{target}-summary.json", summary),
            (f"{target}-bank.parquet", source / output / "oof.parquet"),
            (f"{target}-fits.csv", source / output / "fits.csv"),
        ]:
            inputs[alias] = path
        for component in json.loads(summary.read_text())["config"]["components"]:
            inputs[f"{target}-{component['name']}.parquet"] = source / "inputs" / component["file"]
        targets[target] = {
            "endpoint": target.upper() + "_pIC50_direct_inhibition",
            "bank_oof": f"{target}-bank.parquet",
            "bank_summary": f"{target}-summary.json",
            "bank_fits": f"{target}-fits.csv",
            "method": method,
            "filter": {"method": method},
            "expected_rows": count,
        }
    inputs["cyp2d6-bank.parquet"] = (
        RUNS / "regression/20260926T085959697226Z-0dbf57d4/outputs/oof.parquet"
    )
    inputs["d6-base.parquet"] = (
        RUNS / "regression/20260926T013148480112Z-72cf4682/inputs/control-oof.parquet"
    )
    targets["cyp2d6"] = {
        "endpoint": "CYP2D6_pIC50_direct_inhibition",
        "bank_oof": "cyp2d6-bank.parquet",
        "method": "equal_all",
        "expected_rows": 1457,
        "filter": {
            "split_id": "parent",
            "robust_delta": 0.5,
            "constraint_weight": 0.35,
            "placement": "fold_moment_match",
            "blend_weight": 0.5,
        },
    }
    for item in screens:
        target = item["target"]
        source = RUNS / "regression" / item["run"] / "outputs"
        inputs[f"{target}-trunk.parquet"] = source / "predictions.parquet"
        inputs[f"{target}-cohort.parquet"] = source / "evaluation-cohort.parquet"
        targets[target].update(
            {"predictions": f"{target}-trunk.parquet", "cohort": f"{target}-cohort.parquet"}
        )
    identifier = recorded(
        directory,
        "regression-incumbent-placement",
        "regression",
        {"targets": targets},
        inputs,
        "evaluate_cyp_boltz_incumbent.py",
    )
    (directory / "regression-placement-run.txt").write_text(identifier + "\n")
    if not (directory / "tdi-placement-run.txt").exists():
        subprocess.run(
            ["pixi", "run", "python", "scripts/launch_cyp_boltz_tdi_placement.py", str(directory)],
            cwd=ROOT,
            check=True,
        )
    subprocess.run(
        ["pixi", "run", "python", "scripts/report_cyp_boltz_incumbent.py", str(directory)],
        cwd=ROOT,
        check=True,
    )
    print(f"suite_completed={directory}", flush=True)


if __name__ == "__main__":
    main()
