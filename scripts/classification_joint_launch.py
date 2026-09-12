"""Record all classifier training inputs, including external-only weights."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
base = b / "shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [
    base / "compounds.parquet",
    base / "measurements.parquet",
    b / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    base / "test_mapping.parquet",
    b / "shared/20260912T093252891365Z-271e3020/outputs/observation-targets.parquet",
    b / "shared/20260912T051543210441Z-a3ea2ec3/outputs/best.pt",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Joint TDI Chemprop BCE-only versus auxiliary regression with inner selection and outer refit",
    "--config",
    str(root / "configs/classification-round/joint-v1.json"),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/classification_joint.py"], cwd=root, check=True)
