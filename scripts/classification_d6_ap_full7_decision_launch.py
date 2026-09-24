"""Record full seven-atom AP decision comparisons on the complete UMAP pool."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("variant", choices=["fullfit", "nested"])
variant = parser.parse_args().variant
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
suffix = "decision" if variant == "fullfit" else "nested-umap"
config = root / f"configs/analysis/d6-ap-full7-{suffix}-v1.json"
inputs = {
    "extra-all-pool.parquet": runs / "20260924T033028487589Z-c42e6f27/outputs/all-pool-oof.parquet",
    "ap-full7-all-pool.parquet": runs
    / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
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
    f"Compare complete-pool seven-atom AP {variant} replay with leading atom",
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
