"""Audit retained readouts from the focused-task architecture LGBM screen."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "artifacts/runs/regression/20260916T115927615845Z-1b3e7ef9/outputs"
CONFIG = (
    ROOT
    / "configs/owned-encoder-v2/focused_cyp14_architecture_lgbm_audit_v1.json"
)
files = {
    "w600-d3-ffn3.parquet": "focused_w600_d3_ffn3",
    "w600-d2-ffn3.parquet": "focused_w600_d2_ffn3",
    "w600-d4-ffn3.parquet": "focused_w600_d4_ffn3",
    "w448-d3-ffn3.parquet": "focused_w448_d3_ffn3",
    "w768-d3-ffn3.parquet": "focused_w768_d3_ffn3",
    "w600-d3-ffn2.parquet": "focused_w600_d3_ffn2",
    "confirm-ffn2.parquet": "focused_w600_d3_ffn2_confirm_split",
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose",
    "Audit retained target-specific raw LightGBM readouts of focused CYP14 checkpoints",
    "--config", str(CONFIG),
    "--input", "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
]
for alias, directory in files.items():
    command.extend(
        ["--input", f"{alias}={SCREEN / directory / 'external-oof.parquet'}"]
    )
command.extend(
    [
        "--input",
        "previous-owned.parquet="
        + str(
            ROOT
            / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
            / "external-oof.parquet"
        ),
        "--input",
        "jeremy.parquet="
        + str(
            ROOT
            / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
            / "external-oof.parquet"
        ),
        "--",
        "python",
        "scripts/audit_atomic_readout_grid.py",
    ]
)
subprocess.run(command, cwd=ROOT, check=True)
