"""Launch the recorded raw/PCA CheMeleon TabICL comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw_run = "20260916T060715739725Z-b0ddf3c8"
pca_run = "20260912T071619253798Z-73a6fb7e"
config = root / "configs/regression/chemeleon-raw2048-log2fc-tabicl-audit-v1.json"
inputs = {
    "raw-oof.parquet": root / f"artifacts/runs/regression/{raw_run}/outputs/external-oof.parquet",
    "pca-oof.parquet": root / f"artifacts/runs/regression/{pca_run}/outputs/external-oof.parquet",
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
        "Audit matched raw2048 and PCA256 CheMeleon plus log2fc4 TabICL readouts "
        "and their equal ensemble without incremental-weight search"
    ),
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/audit_regression_chemeleon_raw_tabicl.py"],
    cwd=root,
    check=True,
)
