"""Launch grouped-bootstrap selection audit for block fusion."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
oof = root / "artifacts/runs/classification/20260913T141558951845Z-1bdad261/outputs/oof.parquet"
assert oof.is_file(), oof
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Estimate radius-Chemprop block-fusion selection stability with grouped "
        "bootstrap out-of-bag evaluation"
    ),
    "--config",
    str(root / "configs/classification-round/dual-split-block-fusion-bbc-v1.json"),
    "--input",
    str(oof),
]
subprocess.run(
    command + ["--", "python", "scripts/classification_dual_split_bbc.py"],
    cwd=root,
    check=True,
)
