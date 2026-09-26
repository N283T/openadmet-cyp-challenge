"""Record D6 calibration inventory for previously scored submissions."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
submissions = {
    "s0910.parquet": "20260910T132854718239Z-4766cd9f",
    "s0911.parquet": "20260911T132711430401Z-f3149927",
    "s0912.parquet": "20260912T074802231782Z-eaafa860",
    "s0914a.csv": "20260914T024743176950Z-8ef21e06",
    "s0914b.parquet": "20260914T150145548210Z-1ab27b1d",
    "s0917.csv": "20260917T132706271289Z-193de057",
    "s0918.csv": "20260918T121343664656Z-cf7be0d2",
    "s0924.csv": "20260924T054343781258Z-db4d3cf3",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Align scored D6 submissions and identify affine-only calibration comparisons",
    "--config",
    str(root / "configs/regression/d6-submission-calibration-inventory-v1.json"),
]
for alias, identifier in submissions.items():
    path = root / f"artifacts/submissions/regression/{identifier}/predictions{Path(alias).suffix}"
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_submission_calibration_inventory.py"],
    cwd=root,
    check=True,
)
