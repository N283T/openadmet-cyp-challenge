"""Launch pinned ATTNSOM smoke or full challenge-disjoint training."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["smoke", "full"])
args = parser.parse_args()

root = Path(__file__).resolve().parents[1]
source = Path.home() / "ghq/github.com/dmis-lab/ATTNSOM"
base = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [base / "compounds.parquet"] + [
    source / "cyp_dataset" / f"{isoform}.sdf"
    for isoform in ["1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"]
] + [
    source / filename
    for filename in [
        "main.py",
        "dataset.py",
        "dataset_utils.py",
        "train.py",
        "metrics.py",
        "model.py",
        "LICENSE",
    ]
]
config = root / f"configs/classification-round/attnsom-{args.mode}-v1.json"
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Smoke-test" if args.mode == "smoke" else "Train"
    )
    + " pinned challenge-disjoint ATTNSOM for atom-level SoM representations",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/attnsom_train.py"], cwd=root, check=True)
