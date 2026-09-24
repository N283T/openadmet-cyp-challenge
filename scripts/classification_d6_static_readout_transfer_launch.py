"""Compare the completed D6 LightGBM and ExtraTrees all-pool refits."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/analysis/d6-static-readout-transfer-v1.json"
inputs = {
    "lgbm-all-pool.parquet": root
    / "artifacts/runs/classification/20260924T032014439999Z-c8368f50/outputs/all-pool-oof.parquet",
    "extra-all-pool.parquet": root
    / "artifacts/runs/classification/20260924T032530375340Z-2542e03a/outputs/all-pool-oof.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {*declared["predictions"], declared["assignments"]}:
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
    "Test complete-pool D6 LGBM versus ExtraTrees decision transfer across splits",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_static_full_pool_decision.py"],
    cwd=root,
    check=True,
)
