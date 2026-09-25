"""Launch a recorded matched CV readout of the intermediate CheMeleon states."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
extract_run = "20260925T065503727832Z-1c269e23"
features = root / "artifacts/runs/shared" / extract_run / "outputs"
dataset = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
split = root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
classification = root / "artifacts/runs/classification/20260912T115018743045Z-269e598d/outputs"
sources = {
    **{f"depth-{depth}.parquet": features / f"depth-{depth}.parquet" for depth in (1, 3, 6)},
    "assignments.parquet": split / "assignments.parquet",
    "measurements.parquet": dataset / "measurements.parquet",
    "classification-labels.parquet": classification / "classification-labels.parquet",
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
}
summary = json.loads((features / "summary.json").read_text())
assert summary["final_vs_chemprop_max_abs"] <= 1e-5
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    (
        "Matched UMAP five-fold LightGBM readout of CheMeleon depths 1, 3 and 6 "
        "on four regression and two classification endpoints"
    ),
    "--config",
    str(root / "configs/screen/chemeleon-depth-readout-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_chemeleon_depths.py"])
subprocess.run(command, cwd=root, check=True)
