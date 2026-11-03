"""Launch the frozen MONROE embedding export in the dedicated pixi environment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
checkpoint = Path.home() / "ghq/github.com/blazejba/monroe/checkpoint"
inputs = {
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "monroe-config.json": checkpoint / "config.json",
    "monroe-weights.pt": checkpoint / "weights.pt",
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
        "Export frozen MONROE (blazejba/monroe de075ea) 720-d embeddings for all 6896 "
        "compounds, mean over five seeded ETKDGv3+MMFF94s conformers"
    ),
    "--config",
    str(root / "configs/shared/monroe-embed-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", str(root / ".pixi/envs/monroe/bin/python"), "scripts/monroe_embed.py"],
    cwd=root,
    check=True,
)
