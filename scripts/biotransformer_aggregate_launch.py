"""Launch robust aggregation of the pinned BioTransformer product SDF."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
generation = root / "artifacts/runs/classification/20260914T051230863257Z-abb7661e/outputs"
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    generation / "products.sdf",
]
config = root / "configs/classification-round/biotransformer-aggregate-v1.json"
for path in [config, *inputs]:
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
    "Robustly map and aggregate pinned BioTransformer products without rerunning Java",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/biotransformer_aggregate.py"],
    cwd=root,
    check=True,
)
