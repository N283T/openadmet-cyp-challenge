"""Launch train/test Morgan density and coverage analysis."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
paths = [
    shared / "compounds.parquet",
    shared / "measurements.parquet",
    shared / "test_mapping.parquet",
]
assert all(path.is_file() for path in paths)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Compare hashed-bit positions, radius densities, unhashed environment coverage "
        "and nearest-neighbor support between labeled training and blind test molecules"
    ),
    "--config",
    str(root / "configs/classification-round/fp-train-test-shift-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_fp_train_test_shift.py"],
    cwd=root,
    check=True,
)
