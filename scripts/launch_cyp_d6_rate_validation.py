"""Launch the paired D6 call-count validation."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = {
    "testlike-pool-oof.parquet": runs
    / "20261002T215345688167Z-7920b49a/outputs/oof_scores.parquet",
    "assignments.parquet": runs / "20261002T221918697661Z-2d8d20f5/inputs/assignments.parquet",
    "analog-pairs.parquet": root
    / "artifacts/tracked/d6-analog-expansion-validation-20260924/anchor-analog-pairs.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Paired grouped-bootstrap validation of D6 call counts for the deployed L4 score on "
        "Metric O, Metric M and the pseudo analog-expansion pool"
    ),
    "--config",
    str(root / "configs/classification/d6-rate-validation-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/evaluate_cyp_d6_rate_validation.py"], cwd=root, check=True
)
