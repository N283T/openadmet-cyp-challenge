"""Launch the atomic OOF correlation audit for the reviewed regression banks."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run_inputs = (
    root
    / "artifacts/runs/regression/20260918T020647743819Z-10fdb14b/inputs"
)
config = root / "configs/regression/flat-bank-review-v1.json"
inputs = {path.name: path for path in run_inputs.iterdir() if path.is_file()}
declared = json.loads(config.read_text())
needed = {declared["scorer"]}
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
    "Audit target-specific atomic OOF performance and prediction correlations",
    "--config",
    str(config),
]
for alias in sorted(needed):
    command.extend(["--input", f"{alias}={inputs[alias]}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_bank_correlation_audit.py"],
    cwd=root,
    check=True,
)
