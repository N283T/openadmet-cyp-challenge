"""Launch packaging of the D6-only SMARTCyp live blend candidate."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/classification/20260914T151707304794Z-428271d1/outputs/predictions.parquet",
    root
    / "artifacts/candidates/classification/20260914T020749578109Z-0a0bd4cc/predictions.parquet",
    root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Package audited D6-only 3:1 live/SMARTCyp rank blend and restore current A4; no upload",
    "--config",
    str(root / "configs/classification-round/d6-smartcyp-live-candidate-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_smartcyp_candidate.py"],
    cwd=root,
    check=True,
)
