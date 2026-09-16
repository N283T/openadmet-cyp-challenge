"""Register the audit of selected D6 ranking/kernel/classical blends."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = [
    runs / "20260914T125646427763Z-1707a0e2/outputs/oof.parquet",
    runs / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-ranking-kernel-classical-screen-audit-v1.json"
payload = json.loads(config.read_text())
config.write_text(json.dumps(payload, indent=2) + "\n")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit selected CYP2D6 ranking, kernel, ordinal, spline, and robust blends",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/audit_d6_model_family_screen.py"],
    cwd=root,
    check=True,
)
