"""Launch the compressed-Caruana plus unhashed-radius comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
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
        "Test fixed late-rank additions of the A4 ordinary-plus-unhashed radius-0/1 "
        "block to the two-member compressed Caruana recipe"
    ),
    "--config",
    str(root / "configs/classification-round/caruana-compressed-unhashed-v2.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_caruana_compressed_unhashed.py"],
    cwd=root,
    check=True,
)
