"""Launch the CheMeleon checkpoint similarity analysis (embedding CKA and OOF correlations)."""

import json
import subprocess
from pathlib import Path

from launch_openadmet_depth_screen import CONTROLS, depth_runs

root = Path(__file__).resolve().parents[1]
config = root / "configs/regression/depth-similarity-3a4-v1.json"
depth = json.loads(config.read_text())["depth"]
inputs = {"control-3a4.parquet": CONTROLS["3a4"]}
for name, run in sorted(depth_runs().items()):
    inputs[f"emb-{name}.parquet"] = run / f"outputs/depth-{depth}.parquet"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    f"Similarity of 16 CheMeleon checkpoints at depth {depth}: embedding CKA and OOF correlations",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/regression_depth_similarity.py"], cwd=root, check=True
)
