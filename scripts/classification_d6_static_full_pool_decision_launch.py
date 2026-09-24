"""Record D6 complete-pool call-rate evaluation for two fixed readouts."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/analysis/d6-static-full-pool-decision-v1.json"
inputs = {
    "all-pool-oof.parquet": root
    / "artifacts/runs/classification/20260924T032014439999Z-c8368f50/outputs/all-pool-oof.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {declared["predictions"], declared["assignments"]}:
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
    "Compare full-pool versus labelled-only D6 decisions on UMAP and scaffold refits",
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
