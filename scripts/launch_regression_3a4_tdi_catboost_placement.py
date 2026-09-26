"""Record fixed blends after label-free placement of TDI CatBoost 3A4 OOF."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "scored-bank-oof.parquet": root
    / "artifacts/runs/regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "catboost-oof.parquet": root
    / "artifacts/runs/regression/20260926T132534292013Z-dc82a774/outputs/catboost-oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test label-free fold placement of classification-trained CatBoost peer in scored 3A4 bank",
    "--config",
    str(root / "configs/regression/3a4-tdi-catboost-placement-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_tdi_catboost_placement.py"],
    cwd=root,
    check=True,
)
