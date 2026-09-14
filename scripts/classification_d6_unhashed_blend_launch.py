"""Launch the recorded live-plus-unhashed D6 blend analysis."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs/test.parquet",
    root / "artifacts/runs/classification/20260913T135321291271Z-f79255ae/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T135321291271Z-f79255ae/outputs/test.parquet",
]
for path in inputs:
    assert path.is_file(), path
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Evaluate fixed late-rank blends of live D6 with the dual-split-supported "
        "unhashed radius-2 model, without packaging or upload"
    ),
    "--config",
    str(root / "configs/classification-round/d6-unhashed-blend-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_unhashed_blend.py"],
    cwd=root,
    check=True,
)
