"""Launch the recorded PXR/CYP exact-overlap audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
pxr = Path("/home/nagaet/pxr-iduction-challenge/data")
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    pxr / "train_activity_db.parquet",
    pxr / "phase_1_unblinded_test.parquet",
    pxr / "single_concentration_train.parquet",
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
        "Audit exact standardized-connectivity overlap between all CYP challenge "
        "compounds and prior PXR pEC50 or single-concentration training sources"
    ),
    "--config",
    str(root / "configs/classification-round/pxr-overlap-audit-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_pxr_overlap_audit.py"],
    cwd=root,
    check=True,
)
