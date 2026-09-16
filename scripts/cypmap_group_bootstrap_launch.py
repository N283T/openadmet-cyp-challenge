"""Launch paired group bootstrap for the CYP-MAP mechanism comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
readout = root / "artifacts/runs/classification/20260914T063614433212Z-c0dc3f32/outputs"
inputs = [readout / "oof.parquet"]
config = root / "configs/classification-round/cypmap-group-bootstrap-v1.json"
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Estimate paired grouped intervals for fixed A4 CYP-MAP representation gains",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/cypmap_group_bootstrap.py"], cwd=root, check=True
)
