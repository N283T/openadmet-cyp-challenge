"""Record the fixed incremental CheMeleon depth regression screen."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
sources = {
    "depth-oof.parquet": runs / "20260925T130401763912Z-dcf9b1b8/outputs/oof.parquet",
    "base-1a2.parquet": runs
    / "20260918T034159381060Z-46314b50/outputs/1A2/add_dense_and_fast_pool/oof.parquet",
    "base-2c9.parquet": runs / "20260917T123601175420Z-e1560701/outputs/oof.parquet",
    "base-3a4.parquet": runs / "20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "base-d6-tabpfn.parquet": runs / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "base-d6-lgbm.parquet": runs / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
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
    "Test fixed CheMeleon depth additions to submitted regression OOF recipes",
    "--config",
    str(root / "configs/screen/chemeleon-depth-incremental-regression-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_chemeleon_depth_incremental_regression.py"])
subprocess.run(command, cwd=root, check=True)
