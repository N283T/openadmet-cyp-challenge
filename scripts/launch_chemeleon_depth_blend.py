"""Launch a recorded fixed OOF blend screen of CheMeleon depths."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sources = {
    "oof.parquet": root
    / "artifacts/runs/shared/20260925T065626160287Z-8017d516/outputs/oof.parquet",
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
    "shared",
    "--purpose",
    "Test prespecified equal OOF blends of intermediate and final CheMeleon readouts",
    "--config",
    str(root / "configs/screen/chemeleon-depth-blend-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/chemeleon_depth_blend.py"])
subprocess.run(command, cwd=root, check=True)
