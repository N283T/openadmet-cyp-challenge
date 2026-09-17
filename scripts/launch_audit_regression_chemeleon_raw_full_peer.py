"""Launch the equal-member CheMeleon raw ensemble audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/regression/chemeleon-raw2048-full-peer-audit-v1.json"
run_outputs = root / "artifacts/runs/regression"
inputs = {
    "raw-oof.parquet": run_outputs / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
    "pca-oof.parquet": run_outputs / "20260912T071619253798Z-73a6fb7e/outputs/external-oof.parquet",
    "frozen-oof.parquet": run_outputs
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "scores-oof.parquet": run_outputs
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "finetuned-oof.parquet": run_outputs
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "cyp3a4-top2.parquet": run_outputs
    / "20260912T152832952478Z-3e4de7c8/outputs/predictions.parquet",
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
    (
        "Evaluate raw CheMeleon as a complete equal-weight member of the current "
        "endpoint ensembles without incremental-weight search"
    ),
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/audit_regression_chemeleon_raw_full_peer.py"],
    cwd=root,
    check=True,
)
