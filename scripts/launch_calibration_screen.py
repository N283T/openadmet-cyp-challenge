"""Launch the declared calibration screen with immutable source evidence."""

import subprocess

inputs = [
    "artifacts/runs/regression/20260911T050944039369Z-36d00d9c/outputs/regression-oof-diagnostics.parquet",
    "artifacts/runs/shared/20260910T132141471698Z-343eb5ec/outputs/ensemble_probabilities.parquet",
    "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "artifacts/metric-reference/custom_scoring_functions.py",
    "/home/nagaet/pxr-iduction-challenge/track1_activity/src/importance_weights.py",
    "/home/nagaet/openadmet-pxr-model-report/MODEL_REPORT.md",
    "/home/nagaet/pxr-iduction-challenge/docs/deepresearch/pxr_challenge_regression_report_2026-04-19.md",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Screen 12 PXR-inspired calibration variants using fixed OOF; exploratory, no external labels",
    "--config",
    "configs/calibration/screen-v1.json",
]
for path in inputs:
    cmd += ["--input", path]
subprocess.run(cmd + ["--", "python", "scripts/calibration_screen.py"], check=True)
