"""Launch audit and packaging for the diverse CYP2D6 candidates."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = [
    root / "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs/test.parquet",
    root
    / (
        "artifacts/runs/classification/20260913T103046578961Z-129fbc5d/outputs/"
        "current_d6_caruana_a4_old_count.parquet"
    ),
    root
    / (
        "artifacts/runs/classification/20260912T133350390925Z-effedde5/outputs/"
        "retain_d6_restore_a4.parquet"
    ),
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
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
    "Audit fold stability and package diverse CYP2D6 fixed-rate candidates; no upload",
    "--config",
    str(root / "configs/classification-round/d6-diverse-package-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_diverse_package.py"],
    cwd=root,
    check=True,
)
