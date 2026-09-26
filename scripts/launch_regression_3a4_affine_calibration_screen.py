"""Record affine calibration sensitivity for the scored CYP3A4 bank."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg_path = root / "configs/regression/3a4-affine-calibration-screen-v1.json"
cfg = json.loads(cfg_path.read_text())
inputs = {
    "scored-bank-oof.parquet": root
    / "artifacts/runs/regression"
    / cfg["oof_run"]
    / "outputs/oof.parquet",
    "current.csv": root
    / "artifacts/submissions/regression"
    / cfg["current_submission_id"]
    / "predictions.csv",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Screen 3A4 affine correction on its exactly scored parent bank",
    "--config",
    str(cfg_path),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_affine_calibration_screen.py"],
    cwd=root,
    check=True,
)
