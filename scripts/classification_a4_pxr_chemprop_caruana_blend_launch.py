"""Launch strict PXR Chemprop versus live A4 Caruana blend audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
    root / "artifacts/runs/classification/20260913T153945322549Z-91a3a799/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T153945322549Z-91a3a799/outputs/test.parquet",
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
        "Measure diversity and conservative rank blends between strict PXR Chemprop "
        "A4 models and the deployed six-member Caruana prediction"
    ),
    "--config",
    str(root / "configs/classification-round/a4-pxr-chemprop-caruana-blend-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_pred_bank_caruana_blend.py"],
    cwd=root,
    check=True,
)
