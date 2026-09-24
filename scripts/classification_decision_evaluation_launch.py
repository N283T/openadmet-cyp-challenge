"""Record a CYP2D6 call-rate audit from frozen OOF predictions."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/classification/d6-decision-evaluation-v1.json"
runs = root / "artifacts/runs"
inputs = {
    "candidate-oof.parquet": runs
    / "classification/20260923T171350543926Z-9081bf31/outputs/oof.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
}
declared = json.loads(config.read_text())
if {declared["predictions"], declared["assignments"]} != set(inputs):
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
    "Audit variable D6 call rates for control and AP-mass90 on frozen OOF",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_decision_evaluation.py"],
    cwd=root,
    check=True,
)
