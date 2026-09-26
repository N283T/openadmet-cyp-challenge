"""Record the submitted A4 equal-seven complete-pool call-rate audit."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/classification/a4-equal7-rate-robustness-20260927.json"
oof = root / (
    "artifacts/runs/classification/20260926T150746748801Z-e74748cd/outputs/all-pool-oof.parquet"
)
inputs = {"equal7-full-pool.parquet": oof}
declared = json.loads(config.read_text())
if set(declared["oof"].values()) != set(inputs):
    raise ValueError("Config and input aliases differ")
for path in [config, *inputs.values()]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi", "run", "cyp", "run", "--track", "classification", "--purpose",
    "Select A4 equal-seven call rate outside each full-pool UMAP held fold",
    "--config", str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
completed = subprocess.run(
    command + ["--", "python", "scripts/classification_call_rate_robustness.py"],
    cwd=root, check=False,
)
if completed.returncode:
    raise SystemExit(f"A4 equal-seven rate audit failed with exit code {completed.returncode}")
