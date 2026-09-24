"""Record ExtraTrees versus September 15 control on two complete-pool folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-extra-vs-sep15-full-pool-v1.json"
inputs = {
    "sep15-umap.parquet": runs / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "sep15-scaffold.parquet": runs / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
    "extra.parquet": runs / "20260924T033028487589Z-c42e6f27/outputs/all-pool-oof.parquet",
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
    "Compare ExtraTrees with September 15 control on UMAP and scaffold complete pools",
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
