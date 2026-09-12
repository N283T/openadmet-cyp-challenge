"""Recorded structural/label-route diagnostic with explicit frozen inputs."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
base = b / "shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [
    base / "compounds.parquet",
    base / "measurements.parquet",
    b / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    base / "test_mapping.parquet",
    b / "shared/20260912T093252891365Z-271e3020/outputs/observation-targets.parquet",
    b / "classification/20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "D6 versus3A4 coverage, label-neighbor permutation and assay-route error audit",
    "--config",
    str(root / "configs/classification-round/d6-neighborhood.json"),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_d6_neighborhood.py"], cwd=root, check=True
)
