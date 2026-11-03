"""Launch recorded all-depth exports for the OpenADMET CheMeleon featurizer checkpoints.

Usage: launch_openadmet_chemeleon_depths.py [name ...]; with no names every checkpoint runs, one
recorded run each. The original CheMeleon checkpoint is exported as `chemeleon`.
"""

import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
weights = root / "data/raw/openadmet-chemeleon-featurizers-zenodo-21725564"
original = root / "artifacts/runs/shared/20260910T012630927160Z-e9576045/inputs/001-chemeleon_mp.pt"
checkpoints = {"chemeleon": original}
checkpoints.update({p.name.split("__")[0]: p for p in sorted(weights.glob("*__s42.pt"))})


def launch(name):
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Export all six frozen message-passing depths of the {name} CheMeleon checkpoint",
        "--config",
        str(root / "configs/features/openadmet-chemeleon-depths-v1.json"),
        "--input",
        "compounds.parquet="
        + str(
            root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
        ),
        "--input",
        f"chemeleon_mp.pt={checkpoints[name]}",
        "--",
        "python",
        "scripts/extract_chemeleon_depths.py",
    ]
    subprocess.run(command, cwd=root, check=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or list(checkpoints):
        launch(name)
