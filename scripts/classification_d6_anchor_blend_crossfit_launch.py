"""Launch the recorded D6 anchor blend-weight cross-fit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
blend_oof = (
    root
    / "artifacts/runs/classification/20260918T060729827213Z-8dae5f36/outputs/oof.parquet"
)
if not blend_oof.is_file():
    raise FileNotFoundError(blend_oof)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Cross-fit the fixed CYP2D6 anchor regression blend weights",
    "--config",
    str(root / "configs/classification/d6-anchor-blend-crossfit-v1.json"),
    "--input",
    f"blend-oof.parquet={blend_oof}",
    "--",
    "python",
    "scripts/classification_d6_anchor_blend_crossfit.py",
]
subprocess.run(command, cwd=root, check=True)
