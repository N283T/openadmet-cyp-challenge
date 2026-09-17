"""Record the focused-task checkpoint raw LightGBM screen."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCH = ROOT / "artifacts/runs/shared/20260916T112240065253Z-6543c9c0/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
CONFIG = (
    ROOT
    / "configs/owned-encoder-v2/focused_cyp14_architecture_lgbm_screen_v1.json"
)
inputs = {
    "w600-d3-ffn3.parquet": ARCH / "trial-00-w600_d3_ffn3/external-embeddings.parquet",
    "w600-d2-ffn3.parquet": ARCH / "trial-01-w600_d2_ffn3/external-embeddings.parquet",
    "w600-d4-ffn3.parquet": ARCH / "trial-02-w600_d4_ffn3/external-embeddings.parquet",
    "w448-d3-ffn3.parquet": ARCH / "trial-03-w448_d3_ffn3/external-embeddings.parquet",
    "w768-d3-ffn3.parquet": ARCH / "trial-04-w768_d3_ffn3/external-embeddings.parquet",
    "w600-d3-ffn2.parquet": ARCH / "trial-05-w600_d3_ffn2/external-embeddings.parquet",
    "confirm-ffn2.parquet": ROOT
    / "artifacts/runs/shared/20260916T113640802695Z-f8c7c833/outputs"
    / "external-embeddings.parquet",
    "previous-owned.parquet": ROOT
    / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
    / "external-oof.parquet",
    "jeremy.parquet": ROOT
    / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
    / "external-oof.parquet",
    "measurements.parquet": CHALLENGE / "measurements.parquet",
    "parent-assignments.parquet": ASSIGNMENTS,
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Target-specific raw LightGBM screen of focused CYP14 checkpoints",
    "--config", str(CONFIG),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
command += ["--", "python", "scripts/focused_cyp14_architecture_lgbm_screen.py"]
subprocess.run(command, cwd=ROOT, check=True)
