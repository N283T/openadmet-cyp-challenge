"""Launch the frozen OrbMol embedding export in the dedicated pixi environment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    (
        "Export frozen OrbMol-v2 (orb-models 0.7.0) pooled node features and latent-charge "
        "statistics for all 6896 compounds, mean over five seeded ETKDGv3+MMFF94s conformers"
    ),
    "--config",
    str(root / "configs/shared/orbmol-embed-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", str(root / ".pixi/envs/orbmol/bin/python"), "scripts/orbmol_embed.py"],
    cwd=root,
    check=True,
)
