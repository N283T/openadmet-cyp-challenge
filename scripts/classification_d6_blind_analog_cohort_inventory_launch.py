"""Record blind analog-cohort ID inventory from canonical aliases."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/analysis/d6-blind-analog-cohort-inventory-v1.json"
inputs = {
    "aliases.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet"
}
declared = json.loads(config.read_text())
if set(inputs) != {declared["aliases"]}:
    raise ValueError("Config and input aliases differ")
for path in [config, *inputs.values()]:
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
    "Inventory blind analog-expansion ID block absent from the training cohort",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_blind_analog_cohort_inventory.py"],
    cwd=root,
    check=True,
)
