"""Record uncertainty for the close CYP1A2 whole-member depth candidates."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
sources = {
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
    "parent-assignments.parquet": (
        runs / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "current-ft-lgbm.parquet": (
        runs / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet"
    ),
    "scenario-oof.parquet": runs / "20260925T134831201105Z-04d782a3/outputs/oof.parquet",
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression", "--purpose",
    "Parent-group bootstrap of close CYP1A2 whole-member depth alternatives",
    "--config", str(root / "configs/screen/other-regression-depth-uncertainty-v1.json"),
]
for alias, path in sources.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_d6_depth_atomic_uncertainty.py"])
subprocess.run(command, cwd=root, check=True)
