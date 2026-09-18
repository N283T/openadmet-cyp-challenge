"""Launch the scored-submission-aware v9 placement review."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/regression"
inputs = {
    "old-oof.parquet": run / "20260917T123605451416Z-d9ca2bee/outputs/oof.parquet",
    "new-oof.parquet": run / "20260918T034159381060Z-46314b50/outputs/2D6/v7/oof.parquet",
    "old-raw.csv": run / "20260917T131642179677Z-cf8c2eb0/outputs/candidate.csv",
    "new-raw.csv": run / "20260918T035106389077Z-b97ec048/outputs/candidate-raw.csv",
    "current-placed.csv": run / "20260918T035106389077Z-b97ec048/outputs/candidate-placed.csv",
    "submitted-20260912.parquet": root
    / "artifacts/submissions/regression/20260912T074802231782Z-eaafa860/predictions.parquet",
    "submitted-20260914.csv": root
    / "artifacts/submissions/regression/20260914T024743176950Z-8ef21e06/predictions.csv",
    "submitted-20260917.csv": root
    / "artifacts/submissions/regression/20260917T132706271289Z-193de057/predictions.csv",
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
}
config = root / "configs/regression/submission-aware-placement-v9.json"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Review scored regression placements and package an OOF-guarded D6 delta correction",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/submission_aware_placement_v9.py"],
    cwd=root,
    check=True,
)
