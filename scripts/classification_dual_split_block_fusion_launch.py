"""Launch the recorded dual-split radius/Chemprop block fusion."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/classification/20260913T134843805812Z-0d2682af/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T134843805812Z-0d2682af/outputs/test.parquet",
    root / "artifacts/runs/classification/20260913T135321291271Z-f79255ae/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T135321291271Z-f79255ae/outputs/test.parquet",
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
        "Evaluate prespecified endpoint-specific late fusion of unhashed-radius "
        "and frozen-Chemprop blocks on UMAP and scaffold folds"
    ),
    "--config",
    str(root / "configs/classification-round/dual-split-block-fusion-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_dual_split_block_fusion.py"],
    cwd=root,
    check=True,
)
