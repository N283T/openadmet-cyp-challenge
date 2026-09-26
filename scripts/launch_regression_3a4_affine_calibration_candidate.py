"""Record a 3A4-only correction and its addition to the 1A2/D6 probe."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg_path = root / "configs/regression/3a4-affine-calibration-candidate-v1.json"
cfg = json.loads(cfg_path.read_text())
inputs = {
    "current.csv": root
    / "artifacts/submissions/regression"
    / cfg["source_submission"]
    / "predictions.csv",
    "combined.csv": root
    / "artifacts/submissions/regression"
    / cfg["combined_submission"]
    / "predictions.csv",
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Package conservative 3A4 affine correction alone and with 1A2/D6 probes",
    "--config",
    str(cfg_path),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_affine_calibration_candidate.py"],
    cwd=root,
    check=True,
)
