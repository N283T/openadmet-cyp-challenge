"""Launch the recorded auxiliary-bank versus live Caruana blend audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
    root / "artifacts/runs/classification/20260913T152746349812Z-c510a24e/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T152746349812Z-c510a24e/outputs/test.parquet",
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
        "Measure OOF and blind-rank diversity of the strongest auxiliary prediction "
        "banks against deployed A4 Caruana across conservative blend weights"
    ),
    "--config",
    str(root / "configs/classification-round/a4-pred-bank-caruana-blend-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_pred_bank_caruana_blend.py"],
    cwd=root,
    check=True,
)
