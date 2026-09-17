"""Launch matched PCA-only and raw-for-PCA flat-bank comparisons."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
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
for name in ["flat-ensemble-bank-pca-v1.json", "flat-ensemble-bank-raw-swap-v1.json"]:
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Compare flat regression ensemble methods for matched bank variant {name}",
        "--config",
        str(root / "configs/regression" / name),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"],
        cwd=root,
        check=True,
    )
