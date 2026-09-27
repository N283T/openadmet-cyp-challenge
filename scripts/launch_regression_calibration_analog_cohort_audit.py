"""Record test-selection-like analog sensitivity for 1A2/3A4 calibration."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "1a2-oof.parquet": root
    / "artifacts/runs/regression/20260918T034159381060Z-46314b50"
    / "outputs/1A2/add_dense_and_fast_pool/oof.parquet",
    "3a4-oof.parquet": root
    / "artifacts/runs/regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "analog-pairs.parquet": root
    / "artifacts/tracked/d6-analog-expansion-validation-20260924/anchor-analog-pairs.parquet",
    "parent-assignments.parquet": root
    / "artifacts/runs/regression/20260926T085732408866Z-0167cd1b"
    / "inputs/parent-assignments.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit 1A2 and 3A4 calibration direction on other-CYP selected analog OOF subsets",
    "--config",
    str(root / "configs/regression/calibration-analog-cohort-audit-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_calibration_analog_cohort_audit.py"],
    cwd=root,
    check=True,
)
