"""Launch the recorded A4 live-Caruana versus SoM-representation blend audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/classification/20260912T131834422982Z-d27648b2/outputs/oof.parquet",
    root
    / "artifacts/runs/classification/20260912T131834422982Z-d27648b2/outputs/deployment.parquet",
    root / "artifacts/runs/classification/20260914T051140518895Z-ca43af17/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260914T051140518895Z-ca43af17/outputs/test.parquet",
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
    "Test whether the dual-split SMARTCyp-ATTNSOM A4 block complements live Caruana",
    "--config",
    str(root / "configs/classification-round/a4-som-block-blend-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_block_blend.py"],
    cwd=root,
    check=True,
)
