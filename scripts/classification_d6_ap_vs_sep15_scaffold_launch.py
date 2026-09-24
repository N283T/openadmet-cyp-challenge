"""Record complete-pool scaffold comparison of AP and September 15 control."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("variant", choices=["nested", "rank"])
variant = parser.parse_args().variant
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / f"configs/analysis/d6-ap-vs-sep15-scaffold-{variant}-v1.json"
inputs = {
    "sep15-scaffold.parquet": runs / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
    "ap-scaffold.parquet": runs / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
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
    f"Compare AP and September 15 control on complete scaffold folds: {variant}",
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
