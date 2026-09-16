"""Launch frozen CYP-MAP feature extraction."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["smoke", "full"])
args = parser.parse_args()

root = Path(__file__).resolve().parents[1]
source = Path.home() / "ghq/github.com/qwon135/CYP-MAP"
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
source_files = [
    source / "modules/som_models.py",
    source / "modules/dualgraph/gnn.py",
    source / "modules/dualgraph/conv.py",
    source / "modules/dualgraph/mol.py",
    source / "modules/dualgraph/graph.py",
    source / "modules/ogb/utils/features.py",
]
checkpoints = sorted((source / "ckpt").glob("cypmap_trial*.pt"))
inputs = [shared / "compounds.parquet", *source_files, *checkpoints]
config = root / f"configs/classification-round/cypmap-features-{args.mode}-v1.json"
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)
if len(checkpoints) != 5:
    raise ValueError(f"Expected five CYP-MAP checkpoints, found {len(checkpoints)}")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    f"Extract frozen CYP-MAP reaction-aware features ({args.mode})",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/cypmap_features.py"], cwd=root, check=True)
