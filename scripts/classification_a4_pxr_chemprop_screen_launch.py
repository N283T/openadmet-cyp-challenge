"""Launch the recorded strict-PXR-Chemprop feature screen for A4 TDI."""

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
        "artifacts/runs/classification/20260913T152152397641Z-5177d5d9/outputs/"
        "pxr-classical-features.parquet"
    ),
    root
    / "artifacts/runs/regression/20260911T145924198895Z-373de784/outputs/score-features.parquet",
    root
    / "artifacts/runs/shared/20260912T092919346028Z-bc6e15fc/outputs/substrate-features.parquet",
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
        "Compare strict challenge-disjoint PXR Chemprop predictions and embeddings "
        "alone or with Morgan and classical predicted phenotypes for A4 TDI"
    ),
    "--config",
    str(root / "configs/classification-round/a4-pxr-chemprop-screen-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_pxr_pred_bank.py"],
    cwd=root,
    check=True,
)
