"""Record incremental value of a robust atomic pair to the current 1A2 bank."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "current-oof.parquet": (
        root / "artifacts/tracked/other-regression-depth-members-20260925/ensemble-oof.parquet"
    ),
    "atomic-oof.parquet": (
        root / "artifacts/runs/regression/20260926T015505095663Z-0acd3ef2/outputs/oof.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Check incremental CYP1A2 robust atomic pair against current ten-member OOF",
    "--config", str(root / "configs/regression/1a2-incremental-pair-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_1a2_incremental_pair.py"],
    cwd=root,
    check=True,
)
