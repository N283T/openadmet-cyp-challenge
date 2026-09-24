"""Record the exact September 15/support6/AP-mass90 D6 comparison."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
api = root / "artifacts/api/classification"
config = root / "configs/analysis/d6-three-submission-comparison-v1.json"


def leaderboard_path(run_id: str) -> Path:
    matches = list((api / run_id / "downloads").rglob("*.csv"))
    if len(matches) != 1:
        raise ValueError(f"Expected one CSV in leaderboard snapshot {run_id}")
    return matches[0]


inputs = {
    "control-submission.parquet": runs
    / "20260915T125939829893Z-1d921ffc/outputs/predictions.parquet",
    "support-submission.parquet": runs
    / "20260918T094431957901Z-d535a67b/outputs/predictions.parquet",
    "ap-submission.parquet": runs / "20260923T180205912385Z-1fc93d13/outputs/predictions.parquet",
    "test-diagnostics.parquet": runs
    / "20260918T104416170387Z-9af2247d/outputs/test-diagnostics.parquet",
    "support-oof.parquet": runs / "20260918T085712373599Z-8829d528/outputs/oof.parquet",
    "ap-oof.parquet": runs / "20260923T173655387004Z-97d71dce/outputs/decision-oof.parquet",
    "ap-test.parquet": runs / "20260923T173655387004Z-97d71dce/outputs/decision-test.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "control-leaderboard.csv": leaderboard_path("20260915T130816907483Z-fa8e6053"),
    "support-leaderboard.csv": leaderboard_path("20260918T101216302522Z-4a8b944c"),
    "ap-leaderboard.csv": leaderboard_path("20260923T183846404373Z-b2478be3"),
}
declared = json.loads(config.read_text())
needed = {
    value
    for key, value in declared.items()
    if key
    in [
        "control_submission",
        "support_submission",
        "ap_submission",
        "test_diagnostics",
        "support_oof",
        "ap_oof",
        "ap_test",
        "assignments",
        "control_leaderboard",
        "support_leaderboard",
        "ap_leaderboard",
    ]
}
if set(inputs) != needed:
    raise ValueError("Input aliases differ from config")
for path in [config, *inputs.values()]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Compare why September 15 D6 beat support6 and AP-mass90 on LB",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_three_submission_comparison.py"],
    cwd=root,
    check=True,
)
