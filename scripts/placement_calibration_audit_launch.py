"""Launch independent calibration audit with frozen declared inputs."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
a = b / "regression/20260912T113152922097Z-4a9dd6ab"
n = b / "regression/20260912T113209583937Z-d94acc48"
inputs = [
    a / "inputs/000-all-model-oof.parquet",
    a / "outputs/calibrated-oof.parquet",
    n / "outputs/calibrated-oof.parquet",
    a / "outputs/coefficients.csv",
    a / "outputs/metrics.csv",
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    b / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    n / "inputs/001-custom_scoring_functions.py",
]
cfg = root / "configs/community_reproduction/calibration_audit.json"
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit affine calibration reconstruction and paired D6 group bootstrap500",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/placement_calibration_audit.py"], cwd=root, check=True
)
