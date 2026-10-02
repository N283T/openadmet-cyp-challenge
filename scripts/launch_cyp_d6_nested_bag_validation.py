"""Record the full candidate decision algorithm on untouched outer groups."""

import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/classification/d6-nested-bag-validation-v1.json"
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs"
    shared = runs / "shared/20260909T132734119233Z-a7063ccd/outputs"
    inputs = {
        name: shared / name
        for name in ["compounds.parquet", "measurements.parquet", "aliases.parquet"]
    }
    inputs["assignments.parquet"] = (
        runs / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
    )
    for alias, run in {
        "control-umap.parquet": "20260924T040301085704Z-914bb06c",
        "control-scaffold.parquet": "20260924T041452692950Z-e17778e3",
        "ap-umap.parquet": "20260924T035023873825Z-b719ba35",
        "ap-scaffold.parquet": "20260924T041348966036Z-8db19b4b",
    }.items():
        inputs[alias] = runs / "classification" / run / "outputs/all-pool-oof.parquet"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Validate full five-model D6 decision pipeline before a twelve-hour submission: "
        "untouched outer groups, repeated seeds, aggregation ablations and "
        "historical whole-ensemble diagnostics",
        "--config",
        str(root / args.config),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_zero_start.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
