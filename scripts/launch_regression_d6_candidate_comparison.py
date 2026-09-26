"""Record paired original-versus-robust D6 candidate comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "squared-blend.parquet": runs
    / "regression/20260926T041934073362Z-5e9e8cb8/outputs/oof.parquet",
    "robust-blend.parquet": runs / "regression/20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
    "parent-assignments.parquet": runs
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    "scaffold-assignments.parquet": runs
    / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Pairwise group bootstrap of original and robust D6 candidate OOF",
    "--config",
    str(root / "configs/regression/d6-candidate-comparison-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_candidate_comparison.py"],
    cwd=root,
    check=True,
)
