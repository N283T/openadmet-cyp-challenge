"""Launch the recorded flat CYP3A4 equal-member ensemble audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
config = root / "configs/regression/chemeleon-raw2048-3a4-flat-ensemble-audit-v1.json"
inputs = {
    "caruana-top2.parquet": runs
    / "20260912T152832952478Z-3e4de7c8/outputs/predictions.parquet",
    "frozen-oof.parquet": runs
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "scores-oof.parquet": runs
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "finetuned-oof.parquet": runs
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "pca-oof.parquet": runs
    / "20260912T071619253798Z-73a6fb7e/outputs/external-oof.parquet",
    "raw-oof.parquet": runs
    / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
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
    "Replace nested CYP3A4 pseudo-members with single-layer equal-member recipes",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/audit_regression_3a4_flat_ensemble.py"],
    cwd=root,
    check=True,
)
