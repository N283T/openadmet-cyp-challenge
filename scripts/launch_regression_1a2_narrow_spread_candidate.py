"""Record the 1A2 narrow-spread and combined D6 calibration payloads."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "current.csv": root
    / "artifacts/submissions/regression/20260926T090548323924Z-46188564/predictions.csv",
    "d6-sd90.csv": root
    / "artifacts/submissions/regression/20260926T120124524160Z-01cea54f/predictions.csv",
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
    "Package 1A2-only narrow spread and combined 1A2 plus D6 calibration candidate",
    "--config",
    str(root / "configs/regression/1a2-narrow-spread-candidate-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_1a2_narrow_spread_candidate.py"],
    cwd=root,
    check=True,
)
