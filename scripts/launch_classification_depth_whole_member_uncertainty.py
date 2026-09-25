"""Record labelled TDI whole-member uncertainty."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sources = {
    "scenario-oof.parquet": (
        root / "artifacts/runs/classification/20260925T140800428091Z-fd388a59/outputs/oof.parquet"
    ),
    "umap-assignments.parquet": (
        root / "artifacts/runs/classification/20260912T115018743045Z-269e598d/outputs"
        / "umap-assignments.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "classification", "--purpose",
    "Parent-group bootstrap of selected whole-member depth TDI call changes",
    "--config", str(root / "configs/screen/classification-depth-whole-member-uncertainty-v1.json"),
]
for alias, path in sources.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
command.extend([
    "--", "python", "scripts/evaluate_classification_depth_whole_member_uncertainty.py"
])
subprocess.run(command, cwd=root, check=True)
