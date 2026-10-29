"""Launch the frozen MiniMol embedding export in the dedicated pixi environment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
compounds = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Export frozen MiniMol 1.3.5 (graphium 2.4.7, CPU) 512-d fingerprints for all 6896 compounds",
    "--config",
    str(root / "configs/shared/minimol-embed-v1.json"),
    "--input",
    f"compounds.parquet={compounds}",
]
subprocess.run(
    command + ["--", str(root / ".pixi/envs/minimol/bin/python"), "scripts/minimol_embed.py"],
    cwd=root,
    check=True,
)
