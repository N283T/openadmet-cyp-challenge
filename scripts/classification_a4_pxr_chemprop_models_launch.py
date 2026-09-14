"""Launch fixed strict-PXR-Chemprop embedding model-family screen."""

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
        "Compare classical model families on the fixed challenge-disjoint PXR "
        "Chemprop embedding for A4 TDI under UMAP and scaffold validation"
    ),
    "--config",
    str(root / "configs/classification-round/a4-pxr-chemprop-models-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_pxr_chemprop_models.py"],
    cwd=root,
    check=True,
)
