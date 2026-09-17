"""Launch matched weight-space interpolation of owned joint encoders."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
CONTROL = ROOT / "artifacts/runs/regression/20260916T151717448226Z-ea34957f/outputs"
FOCUSED = ROOT / "artifacts/runs/regression/20260916T150348818193Z-aa2923e8/outputs"
CONFIG = ROOT / "configs/external-training/joint-encoder-weight-interpolation-v1.json"

inputs = {
    "compounds.parquet": CHALLENGE / "compounds.parquet",
    "measurements.parquet": CHALLENGE / "measurements.parquet",
    "parent-assignments.parquet": ASSIGNMENTS,
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
    "continuous-base.pt": CONTROL
    / "trial-00-joint_continuous_focus1_control/checkpoint.pt",
    "cyp-base.pt": CONTROL / "trial-01-joint_cyp_focus1_control/checkpoint.pt",
    "all-base.pt": CONTROL / "trial-02-joint_all_focus1_control/checkpoint.pt",
    "continuous-1a2-focus2.pt": FOCUSED
    / "trial-00-joint_continuous_1a2_focus2/checkpoint.pt",
    "continuous-2c9-focus2.pt": FOCUSED
    / "trial-02-joint_continuous_2c9_focus2/checkpoint.pt",
    "cyp-2d6-focus2.pt": FOCUSED / "trial-04-joint_cyp_2d6_focus2/checkpoint.pt",
    "all-3a4-focus2.pt": FOCUSED / "trial-06-joint_all_3a4_focus2/checkpoint.pt",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test single-encoder weight interpolation between matched joint models",
    "--config",
    str(CONFIG),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/interpolate_joint_encoder_weights.py"],
    cwd=ROOT,
    check=True,
)
