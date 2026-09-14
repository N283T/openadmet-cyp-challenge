"""Launch cross-fit and grouped-bootstrap validation of auxiliary A4 blends."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/classification/20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    root / "artifacts/runs/classification/20260913T152746349812Z-c510a24e/outputs/oof.parquet",
    root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
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
        "Validate conservative auxiliary A4 Caruana blends by leave-one-fold-out "
        "weight selection and UMAP-group bootstrap intervals"
    ),
    "--config",
    str(root / "configs/classification-round/a4-pred-bank-blend-validation-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_pred_bank_blend_validation.py"],
    cwd=root,
    check=True,
)
