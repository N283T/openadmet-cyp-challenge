"""Record AP versus September 15 control on complete UMAP folds."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("variant", choices=["fullfit", "nested", "rank-sensitivity"])
variant = parser.parse_args().variant
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / f"configs/analysis/d6-ap-vs-sep15-{variant}-full-pool-v1.json"
if variant == "rank-sensitivity":
    config = root / "configs/analysis/d6-ap-vs-sep15-rank-sensitivity-v1.json"
inputs = {
    "sep15-all-pool.parquet": runs / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "ap-all-pool.parquet": runs / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
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
    f"Compare complete-pool AP {variant} with September 15 D6 control analog",
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
