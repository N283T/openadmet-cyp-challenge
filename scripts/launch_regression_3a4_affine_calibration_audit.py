"""Record grouped uncertainty for candidate 3A4 affine corrections."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg_path = root / "configs/regression/3a4-affine-calibration-audit-v1.json"
cfg = json.loads(cfg_path.read_text())
inputs = {
    "scored-bank-oof.parquet": root
    / "artifacts/runs/regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "parent-assignments.parquet": root
    / "artifacts/runs/regression/20260926T085732408866Z-0167cd1b/inputs/parent-assignments.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit series-level uncertainty for scored 3A4 affine calibration candidates",
    "--config",
    str(cfg_path),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_affine_calibration_audit.py"],
    cwd=root,
    check=True,
)
