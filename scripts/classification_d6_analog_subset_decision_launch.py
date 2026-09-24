"""Record D6 decisions within the pseudo analog-expansion validation pool."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-analog-subset-decision-v1.json"
inputs = {
    "analog-pairs.parquet": runs
    / "20260924T044113708854Z-1294f049/outputs/anchor-analog-pairs.parquet",
    "assignments.parquet": shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "sep15-umap.parquet": runs / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "ap-umap.parquet": runs / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
    "sep15-scaffold.parquet": runs / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
    "ap-scaffold.parquet": runs / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
}
declared = json.loads(config.read_text())
expected = {
    declared["analogs"],
    declared["assignments"],
    *declared["predictions"]["umap"],
    *declared["predictions"]["scaffold"],
}
if set(inputs) != expected:
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
    "Compare September 15 and AP within a pseudo analog test pool",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_analog_subset_decision.py"],
    cwd=root,
    check=True,
)
