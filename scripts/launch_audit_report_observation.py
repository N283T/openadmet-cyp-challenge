"""Launch a declared-input independent audit of the observation-model experiment."""

import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    run = "20260912T093601580048Z-ec2b7210"
    output = root / f"artifacts/runs/regression/{run}/outputs"
    config = root / "configs/report-observation/independent-audit-v3.json"
    with config.open("x") as stream:
        json.dump(
            dict(model_run=run, arms=["independent", "hill_fixed", "hill_learned"]),
            stream,
            indent=2,
        )
    inputs = [
        output / name
        for name in [
            "observation-oof.parquet",
            "memberships.parquet",
            "training-history.parquet",
            "summary.json",
            "metrics.csv",
        ]
    ]
    inputs += sorted(output.glob("*.pt")) + sorted(output.glob("*-embeddings.parquet"))
    inputs += [
        root
        / "artifacts/runs/regression/20260912T083627802028Z-ae562a59/outputs"
        / "focus_2d6_combined_50pct.parquet",
        root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet",
        root
        / "artifacts/runs/shared/20260912T093252891365Z-271e3020/outputs"
        / "observation-targets.parquet",
        root
        / "artifacts/runs/shared/20260912T055819206227Z-c060d197/outputs"
        / "external-embeddings.parquet",
        root
        / "artifacts/runs/shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ]
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Independent observation-model audit: official OOF, splits, train-only scaling, "
        "epoch selection and CPU reconstruction",
        "--config",
        str(config),
    ]
    for path in inputs:
        assert path.is_file(), path
        command += ["--input", str(path)]
    result = subprocess.run(
        command + ["--", "python", "scripts/audit_report_observation.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    with (root / "artifacts/launches/audit_report_observation_v3.log").open("x") as stream:
        stream.write(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()


if __name__ == "__main__":
    main()
