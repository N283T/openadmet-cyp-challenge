"""Record AP six-atom versus top-atom complete-pool decision comparisons."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("variant", choices=["fullfit", "nested"])
variant = parser.parse_args().variant
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = (
    root
    / f"configs/analysis/d6-ap-top6-{'decision' if variant == 'fullfit' else 'nested-umap'}-v1.json"
)
inputs = {
    "extra-all-pool.parquet": runs / "20260924T033028487589Z-c42e6f27/outputs/all-pool-oof.parquet",
    "ap-top6-all-pool.parquet": runs
    / "20260924T034601208462Z-2792c841/outputs/all-pool-oof.parquet",
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
    f"Compare AP six-atom {variant} proxy with leading atom on complete pools",
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
