"""Launch grouped-bootstrap selection audit for the radius comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
oof = root / "artifacts/runs/classification/20260913T135321291271Z-f79255ae/outputs/oof.parquet"
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
        "Estimate endpoint-specific radius candidate selection stability with "
        "grouped bootstrap out-of-bag evaluation"
    ),
    "--config",
    str(root / "configs/classification-round/dual-split-radius-bbc-v1.json"),
    "--input",
    str(oof),
]
subprocess.run(
    command + ["--", "python", "scripts/classification_dual_split_bbc.py"],
    cwd=root,
    check=True,
)
