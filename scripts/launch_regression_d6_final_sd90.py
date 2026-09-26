"""Record the conservative final D6 SD-0.90 alternative."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "candidate.csv": root
    / "artifacts/runs/regression/20260926T090408836539Z-239369e6/outputs/candidate.csv",
    "current.csv": root
    / "artifacts/submissions/regression/20260924T054343781258Z-db4d3cf3/predictions.csv",
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
    "Prepare D6 final-spread restoration at previously scored moments",
    "--config",
    str(root / "configs/regression/d6-final-sd90-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_final_sd90.py"],
    cwd=root,
    check=True,
)
