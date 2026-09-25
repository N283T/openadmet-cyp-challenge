"""Re-evaluate frozen depth views on the parent-retained regression CV."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
features = root / "artifacts/runs/shared/20260925T065503727832Z-1c269e23/outputs"
sources = {
    **{f"depth-{depth}.parquet": features / f"depth-{depth}.parquet" for depth in (1, 3, 6)},
    "parent-assignments.parquet": (
        root
        / "artifacts/runs/regression/20260912T073049029140Z-0c26dd06/inputs"
        / "005-parent-assignments.parquet"
    ),
    "measurements.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Re-evaluate CheMeleon depths on the current parent-retained regression split",
    "--config",
    str(root / "configs/screen/chemeleon-depth-parent-readout-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_chemeleon_depths.py"])
subprocess.run(command, cwd=root, check=True)
