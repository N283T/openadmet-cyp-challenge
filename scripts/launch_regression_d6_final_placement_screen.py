"""Record D6 final-placement comparisons for the robust-bound candidate."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
path = root / "artifacts/runs/regression/20260926T085959697226Z-0dbf57d4/outputs/oof.parquet"
if not path.is_file():
    raise FileNotFoundError(path)
subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Test label-free final placement preserving scored D6 distribution",
        "--config",
        str(root / "configs/regression/d6-final-placement-screen-v1.json"),
        "--input",
        f"blend-oof.parquet={path}",
        "--",
        "python",
        "scripts/regression_d6_final_placement_screen.py",
    ],
    cwd=root,
    check=True,
)
