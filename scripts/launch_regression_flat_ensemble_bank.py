"""Launch the first regression atomic-bank ensemble comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
config = root / "configs/regression/flat-ensemble-bank-v1.json"
inputs = {
    "jeremy-frozen-oof.parquet": runs
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "jeremy-scores-oof.parquet": runs
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "jeremy-finetuned-oof.parquet": runs
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "chemeleon-pca-oof.parquet": runs
    / "20260912T071619253798Z-73a6fb7e/outputs/external-oof.parquet",
    "chemeleon-raw-oof.parquet": runs
    / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
    "current-control.parquet": runs / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
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
    "Compare auditable flat ensemble methods on an atomic regression prediction bank",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"], cwd=root, check=True
)
