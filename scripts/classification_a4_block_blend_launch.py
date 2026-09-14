"""Launch the recorded A4 Caruana/block blend comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/classification/20260912T131834422982Z-d27648b2/outputs/oof.parquet",
    root
    / "artifacts/runs/classification/20260912T131834422982Z-d27648b2/outputs/deployment.parquet",
    root / "artifacts/runs/classification/20260913T141558951845Z-1bdad261/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T141558951845Z-1bdad261/outputs/test.parquet",
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
        "Test whether the dual-split-supported A4 radius-Chemprop block "
        "complements the live Caruana rank at 249 calls"
    ),
    "--config",
    str(root / "configs/classification-round/a4-block-blend-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_block_blend.py"],
    cwd=root,
    check=True,
)
