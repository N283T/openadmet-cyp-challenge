"""Record conservative control/AP rank blend on UMAP and scaffold folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/classification/d6-control-ap25-full-pool-v1.json"
inputs = {
    "sep15-umap.parquet": runs / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "ap-umap.parquet": runs / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
    "sep15-scaffold.parquet": runs / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
    "ap-scaffold.parquet": runs / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {*declared["predictions"]["umap"], *declared["predictions"]["scaffold"]}:
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
    "Compose 75% September 15 plus 25% AP ranks on complete UMAP and scaffold folds",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_control_ap25_full_pool.py"],
    cwd=root,
    check=True,
)
