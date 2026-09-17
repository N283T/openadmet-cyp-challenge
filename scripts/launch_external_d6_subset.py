"""Launch a recorded direct-CYP2D6 external observation subset."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "artifacts/runs/shared/20260916T072819161332Z-a374a545/outputs"
CONFIG = ROOT / "configs/external-training/owned-encoder-d6-focused-v1.json"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Prepare direct-CYP2D6-only external supervision for a focused project-owned encoder",
    "--config",
    str(CONFIG),
    "--input",
    f"external-compounds.parquet={BASE / 'external-compounds.parquet'}",
    "--input",
    f"observations.parquet={BASE / 'observations.parquet'}",
    "--",
    "python",
    "scripts/external_observation_subset.py",
]
subprocess.run(command, cwd=ROOT, check=True)
