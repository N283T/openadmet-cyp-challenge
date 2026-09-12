"""Record feature assembly from every original Mordred shard and shared compounds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/external-tree/fp_mordred_baseline_v1.json"
with config.open("x") as stream:
    json.dump(
        dict(
            dataset_run="20260909T132734119233Z-a7063ccd",
            mordred_run="20260909T132918829974Z-bc81af58",
            purpose="Raw feature assembly for baseline-equivalent futureCaruana arms",
        ),
        stream,
        indent=2,
    )
files = [root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"]
files += sorted(
    (root / "artifacts/runs/shared/20260909T132918829974Z-bc81af58/outputs").glob(
        "mordred-*.parquet"
    )
)
assert len(files) > 1 and all(path.is_file() for path in files)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Assemble baseline Morgan2048 plus original Mordred; no model or preprocessing fit",
    "--config",
    str(config),
]
for path in files:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/prepare_fp_mordred.py"], cwd=root, check=True)
