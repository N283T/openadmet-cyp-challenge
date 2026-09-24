"""Record hypothetical D6 shift-interval scoring-mask sensitivity."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-ap-scoring-mask-sensitivity-v1.json"
inputs = {
    "decisions.parquet": runs
    / "20260924T042633435522Z-c55d672a/outputs/decisions-by-compound.parquet",
    "assay-route.parquet": runs / "20260912T125901882388Z-4d006374/outputs/assay-route.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {declared["decisions"], declared["assay_route"], declared["assignments"]}:
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
    "Test AP versus September 15 under a hypothetical shift-interval scoring mask",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_scoring_mask_sensitivity.py"],
    cwd=root,
    check=True,
)
