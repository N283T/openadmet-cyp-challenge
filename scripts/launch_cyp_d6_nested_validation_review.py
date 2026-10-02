"""Record interpretation diagnostics of a completed nested validation run."""

import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_id")
    parser.add_argument("--comparison-run-id")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    outputs = root / "artifacts/runs/classification" / args.run_id / "outputs"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Review actual nested D6 decisions, seed stability and label-regime errors "
        "before judging the value of a twelve-hour submission opportunity",
        "--config",
        str(root / "configs/classification/d6-nested-validation-review-v1.json"),
    ]
    for name in [
        "nested-oof.parquet",
        "source-audit.parquet",
        "threshold-scale-audit.csv",
        "leakage-and-operating-points.csv",
        "summary.json",
        "aggregate-metrics.csv",
        "whole-ensemble-ranking-diagnostics.csv",
        "paired-group-intervals.csv",
    ]:
        path = outputs / name
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{name}={path}"])
    if args.comparison_run_id:
        path = (
            root
            / "artifacts/runs/classification"
            / args.comparison_run_id
            / "outputs/nested-oof.parquet"
        )
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"comparison-oof.parquet={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/review_cyp_d6_nested_validation.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
