"""Launch fixed-rate comparison of old and PXR-extended nested Caruana."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
old = root / "artifacts/runs/classification/20260912T141520708830Z-97df4c71/outputs"
new = root / "artifacts/runs/classification/20260913T161419718564Z-2e6c1f38/outputs"
inputs = [
    old / "oof.parquet",
    old / "deployment.parquet",
    new / "oof.parquet",
    new / "deployment.parquet",
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
        "Compare original and strict-PXR-extended nested Caruana recipes at the "
        "same foldwise 33.2 percent OOF and top-249 blind call rates"
    ),
    "--config",
    str(root / "configs/classification-round/caruana-pxr-fixed-rate-audit-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_caruana_pxr_fixed_rate_audit.py"],
    cwd=root,
    check=True,
)
