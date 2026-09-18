"""Launch the target-specific duplicate pruning and bank swap review."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs_dir = runs / "20260918T015803413341Z-7286be9d/inputs"
screen = runs / "20260918T012518053809Z-5c46e72d/outputs"
config = root / "configs/regression/flat-bank-review-v1.json"
inputs = {path.name: path for path in inputs_dir.iterdir() if path.is_file()}
inputs.update(
    {
        "lgbm-jeremy-all-log2fc.parquet": screen / "jeremy_all_log2fc-oof.parquet",
        "lgbm-jeremy-anchor-all-log2fc.parquet": screen
        / "jeremy_anchor_all_log2fc-oof.parquet",
    }
)
declared = json.loads(config.read_text())
needed = {declared["scorer"], declared["current_control"]}
for pool in declared["pools"].values():
    needed.update(file for _, file in pool["components"])
missing = needed - set(inputs)
if missing:
    raise FileNotFoundError(sorted(missing))

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Review target-specific flat banks by prediction redundancy, pruning and swaps",
    "--config",
    str(config),
]
for alias in sorted(needed):
    command.extend(["--input", f"{alias}={inputs[alias]}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_flat_bank_review.py"],
    cwd=root,
    check=True,
)
