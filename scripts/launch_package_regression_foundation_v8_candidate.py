"""Package and register the regression foundation-v8 submission candidate."""

import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
reg = root / "artifacts/runs/regression"
shared = root / "artifacts/runs/shared"
inputs = {
    "aliases.parquet": shared / "20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "previous-submission.csv": reg / "20260917T131642179677Z-cf8c2eb0/outputs/candidate-placed.csv",
    "previous-atomic.parquet": reg
    / "20260917T131642179677Z-cf8c2eb0/outputs/atomic-test-predictions.parquet",
    "new-atomic.parquet": reg / "20260918T034555328092Z-a105b341/outputs/test-predictions.parquet",
    "full-chemeleon-2c9.parquet": reg
    / "20260918T034905633149Z-21ad1846/outputs/test-predictions.parquet",
    "2c9-fits.csv": reg / "20260918T034159381060Z-46314b50/outputs/2C9/add_dense_icl/fits.csv",
    "2d6-fits.csv": reg / "20260918T034159381060Z-46314b50/outputs/2D6/v7/fits.csv",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Package all-target foundation-v8 flat-bank regression submission candidate",
    "--config",
    str(root / "configs/regression/foundation-v8-submission-candidate.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
result = subprocess.run(
    command + ["--", "python", "scripts/package_regression_foundation_v8_candidate.py"],
    cwd=root,
    text=True,
    capture_output=True,
)
print(result.stdout, flush=True)
print(result.stderr, flush=True)
result.check_returncode()
run_ids = re.findall(r"run_id=(\S+)", result.stdout)
if len(run_ids) != 1:
    raise ValueError(run_ids)
run_id = run_ids[0]
subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "candidate",
        "--track",
        "regression",
        "--file",
        str(reg / run_id / "outputs/candidate-placed.csv"),
        "--reference",
        str(root / "data/raw/cyp-challenge-TEST-BLINDED.csv"),
        "--run",
        run_id,
        "--notes",
        (
            "Foundation-v8 target-specific flat banks. All four targets change; "
            "placement moments match the prior submitted candidate. No test labels used."
        ),
    ],
    cwd=root,
    check=True,
)
