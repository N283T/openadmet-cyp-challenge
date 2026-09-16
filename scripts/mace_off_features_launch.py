"""Launch pinned MACE-OFF23 embedding extraction."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["smoke", "full"])
args = parser.parse_args()

root = Path(__file__).resolve().parents[1]
model = (
    Path.home() / "ghq/github.com/ACEsuit/mace-off/mace_off23/MACE-OFF23_small.model"
)
compounds = (
    root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
)
config = root / f"configs/classification-round/mace-off-{args.mode}-v1.json"
for path in (model, compounds, config):
    if not path.is_file():
        raise FileNotFoundError(path)

subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        ("Smoke-test" if args.mode == "smoke" else "Extract")
        + " pinned ligand-only MACE-OFF23 embeddings for TDI classification",
        "--config",
        str(config),
        "--input",
        str(model),
        "--input",
        str(compounds),
        "--",
        "python",
        "scripts/mace_off_features.py",
    ],
    cwd=root,
    check=True,
)
