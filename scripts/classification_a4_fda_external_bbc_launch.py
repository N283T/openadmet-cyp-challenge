"""Launch grouped-bootstrap selection audit for the A4 FDA prediction block."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
oof = root / "artifacts/runs/classification/20260913T143341518156Z-5463a4dc/outputs/oof.parquet"
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
        "Estimate grouped-bootstrap out-of-bag selection stability for the "
        "independent FDA CYP3A4 TDI prediction-block family"
    ),
    "--config",
    str(root / "configs/classification-round/a4-fda-external-bbc-v1.json"),
    "--input",
    str(oof),
    "--",
    "python",
    "scripts/classification_dual_split_bbc.py",
]
subprocess.run(command, cwd=root, check=True)
