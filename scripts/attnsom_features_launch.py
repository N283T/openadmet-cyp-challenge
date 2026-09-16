"""Launch challenge-wide ATTNSOM SoM and embedding extraction."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = Path.home() / "ghq/github.com/dmis-lab/ATTNSOM"
train = root / "artifacts/runs/classification/20260914T044928216258Z-221787f2/outputs"
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [
    shared / "compounds.parquet",
    source / "dataset_utils.py",
    source / "model.py",
    train / "checkpoints/config.json",
    *[train / f"checkpoints/fold{fold}_best.pt" for fold in range(1, 6)],
]
config = root / "configs/classification-round/attnsom-features-v1.json"
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
    "Extract challenge-disjoint ATTNSOM atom-level SoM scores and molecular embeddings",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/attnsom_features.py"], cwd=root, check=True)
