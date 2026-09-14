"""Launch the recorded A4 Caruana subset-distillation audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
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
        "Distill the deployed six-member A4 Caruana prediction into every smaller "
        "member subset using OOF-only simplex weights, then measure blind ranking "
        "and top-249 call fidelity without blind labels"
    ),
    "--config",
    str(root / "configs/classification-round/caruana-compression-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_caruana_compression.py"],
    cwd=root,
    check=True,
)
