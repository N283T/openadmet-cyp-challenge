"""Launch nested strict-PXR-Chemprop predictions for Caruana selection."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    root
    / (
        "artifacts/runs/classification/20260913T153545798308Z-519bbb46/outputs/"
        "pxr-chemprop-features.parquet"
    ),
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
        "Generate matched outer, rotated-inner, and fold-test probabilities for "
        "strict PXR Chemprop readouts so Caruana can select them without OOF leakage"
    ),
    "--config",
    str(root / "configs/classification-round/pxr-chemprop-nested-bank-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_pxr_chemprop_nested_bank.py"],
    cwd=root,
    check=True,
)
