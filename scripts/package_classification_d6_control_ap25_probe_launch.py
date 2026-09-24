"""Record the conservative D6/AP blind probe package without uploading."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
config = root / "configs/submission-classification/d6-control-ap25-probe-v1.json"
inputs = {
    "base-submission.parquet": runs
    / "classification/20260915T125939829893Z-1d921ffc/outputs/predictions.parquet",
    "decision-test.parquet": runs
    / "classification/20260923T173655387004Z-97d71dce/outputs/decision-test.parquet",
    "aliases.parquet": runs / "shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
}
declared = json.loads(config.read_text())
if set(inputs) != {
    declared["base_submission"],
    declared["decision_test"],
    declared["aliases"],
    declared["reference"],
}:
    raise ValueError("Config and input aliases differ")
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
    "Prepare unsubmitted four-swap conservative D6 control/AP leaderboard probe",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/package_classification_d6_control_ap25_probe.py"],
    cwd=root,
    check=True,
)
