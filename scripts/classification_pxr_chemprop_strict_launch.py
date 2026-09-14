"""Launch challenge-disjoint PXR Chemprop training and feature extraction."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    Path("/home/nagaet/pxr-iduction-challenge/data/single_concentration_train.parquet"),
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
        "Retrain a three-seed PXR log2fc Chemprop ensemble after excluding every "
        "current CYP connectivity, then extract safe predictions and embeddings"
    ),
    "--config",
    str(root / "configs/classification-round/pxr-chemprop-strict-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_pxr_chemprop_strict.py"],
    cwd=root,
    check=True,
)
