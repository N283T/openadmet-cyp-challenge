"""Recorded matched-pool follow-up to neighborhood diagnostic."""

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
    b / "classification/20260912T125901882388Z-4d006374/outputs/assay-route.parquet",
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
    "Control neighbor-pool size and decompose selected ensemble assay-regime confusion",
    "--config",
    str(root / "configs/classification-round/d6-neighborhood.json"),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_d6_coverage_control.py"], cwd=root, check=True
)
