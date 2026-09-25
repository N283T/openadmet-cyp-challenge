"""Record intermediate extraction from the current FT600 ensemble member."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ft = root / "artifacts/runs/regression/20260912T060538153709Z-8d434eee/outputs"
sources = {
    "compounds.parquet": (
        root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs" / "compounds.parquet"
    ),
}
for fold in range(5):
    sources[f"fold{fold}-best.pt"] = ft / f"fold{fold}/best.pt"
    sources[f"fold{fold}-final.parquet"] = ft / f"fold{fold}/external-embeddings.parquet"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Extract fold-safe Jeremy FT600 message-passing depth embeddings",
    "--config",
    str(root / "configs/screen/jeremy-ft-depth-extraction-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/extract_jeremy_ft_depths.py"])
subprocess.run(command, cwd=root, check=True)
