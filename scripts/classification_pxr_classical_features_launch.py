"""Launch the recorded, challenge-disjoint classical PXR feature bank."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
pxr = Path("/home/nagaet/pxr-iduction-challenge/data")
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    pxr / "train_activity_db.parquet",
    pxr / "phase_1_unblinded_test.parquet",
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
        "Train diverse classical PXR pEC50 models after excluding every current CYP "
        "connectivity, validate them by scaffold groups, and predict all CYP compounds"
    ),
    "--config",
    str(root / "configs/classification-round/pxr-classical-feature-bank-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_pxr_classical_features.py"],
    cwd=root,
    check=True,
)
