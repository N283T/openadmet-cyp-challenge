"""Launch the recorded FDA CYP3A4 TDI supplement fetch."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Archive and inventory the primary FDA CYP3A4 TDI supplementary dataset "
        "without deriving challenge labels"
    ),
    "--config",
    str(root / "configs/external/fda-cyp3a4-tdi-source-v1.json"),
    "--",
    "python",
    "scripts/fetch_fda_tdi_supplement.py",
]
subprocess.run(command, cwd=root, check=True)
