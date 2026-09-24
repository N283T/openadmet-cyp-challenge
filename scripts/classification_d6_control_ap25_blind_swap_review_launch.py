"""Record blind swap review for the conservative unsubmitted D6 probe."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-control-ap25-blind-swap-review-v1.json"
inputs = {
    "control.parquet": runs / "20260915T125939829893Z-1d921ffc/outputs/predictions.parquet",
    "candidate.parquet": runs / "20260924T043652241766Z-85d10ccf/outputs/predictions.parquet",
    "decision.parquet": runs / "20260924T043652241766Z-85d10ccf/outputs/decision-test.parquet",
    "diagnostics.parquet": runs / "20260923T190811462293Z-4f4a92f2/outputs/threeway-test.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {declared[k] for k in ["control", "candidate", "decision", "diagnostics"]}:
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
    "Review four added and four removed blind D6 calls in unsubmitted 75/25 probe",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_control_ap25_blind_swap_review.py"],
    cwd=root,
    check=True,
)
