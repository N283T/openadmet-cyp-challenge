"""Launch three-way strict PXR Chemprop readout blend screen."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/classification/20260913T153945322549Z-91a3a799/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T153945322549Z-91a3a799/outputs/test.parquet",
    root / "artifacts/runs/classification/20260913T154704125936Z-8e8402cb/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T154704125936Z-8e8402cb/outputs/test.parquet",
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
        "Test three-way rank blends of Morgan plus strict PXR Chemprop LightGBM "
        "with complementary ExtraTrees and RBF-SVM readouts"
    ),
    "--config",
    str(root / "configs/classification-round/a4-pxr-chemprop-model-blend-v2.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_pxr_chemprop_model_blend.py"],
    cwd=root,
    check=True,
)
