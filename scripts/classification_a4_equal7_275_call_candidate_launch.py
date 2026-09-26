"""Record a prospective count-only A4 leaderboard candidate."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification"
config = root / "configs/classification/a4-equal7-275-call-candidate-20260927.json"
inputs = {
    "current-candidate.parquet": root
    / "artifacts/candidates/classification/20260915T130001184144Z-2248c77c/predictions.parquet",
    "submitted-test-scores.parquet": run
    / "20260915T125939829893Z-1d921ffc/outputs/test-scores.parquet",
    "test-aliases.parquet": run
    / "20260915T125939829893Z-1d921ffc/inputs/006-aliases.parquet",
    "blind-reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "a4-rate-curve.csv": run / "20260926T150803596078Z-4d566758/outputs/rate-curve.csv",
}
declared = json.loads(config.read_text())
if set(inputs) != {
    declared[key]
    for key in ("current_candidate", "submitted_scores", "aliases", "reference", "rate_curve")
}:
    raise ValueError("Config and input aliases differ")
for path in [config, *inputs.values()]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi", "run", "cyp", "run", "--track", "classification", "--purpose",
    "Prepare prospective A4 equal-seven top-275 count-only candidate with D6 unchanged",
    "--config", str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
completed = subprocess.run(
    command + ["--", "python", "scripts/classification_a4_equal7_275_call_candidate.py"],
    cwd=root, check=False,
)
if completed.returncode:
    raise SystemExit(f"A4 count-only candidate failed with exit code {completed.returncode}")
