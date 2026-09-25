"""Launch a recorded CheMeleon intermediate-depth export."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Compare frozen CheMeleon intermediate message-passing depths with the current final embedding",
    "--config",
    str(root / "configs/features/chemeleon-depth-probe-v1.json"),
    "--input",
    "compounds.parquet="
    + str(
        root
        / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
    ),
    "--input",
    "chemeleon_mp.pt="
    + str(
        root
        / "artifacts/runs/shared/20260910T012630927160Z-e9576045/inputs/001-chemeleon_mp.pt"
    ),
    "--",
    "python",
    "scripts/extract_chemeleon_depths.py",
]
subprocess.run(command, cwd=root, check=True)
