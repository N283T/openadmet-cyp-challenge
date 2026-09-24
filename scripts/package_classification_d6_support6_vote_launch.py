"""Launch packaging for the user-approved D6 support6 leaderboard probe."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "base-submission.parquet": runs
    / "classification/20260915T125939829893Z-1d921ffc/outputs/predictions.parquet",
    "consensus-test.parquet": runs
    / "classification/20260918T085712373599Z-8829d528/outputs/test.parquet",
    "aliases.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
}
for path in inputs.values():
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
    "Package user-approved D6 support6 vote with the scored A4 vector unchanged",
    "--config",
    str(
        root
        / "configs/submission-classification/d6-support6-vote-a4-unchanged-v1.json"
    ),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/package_classification_d6_support6_vote.py"],
    cwd=root,
    check=True,
)
