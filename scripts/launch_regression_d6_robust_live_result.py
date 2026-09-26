"""Record pre/post Live D6 and macro metrics for the submitted robust blend."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
api = root / "artifacts/api/regression"


def single_csv(run: str) -> Path:
    found = list((api / run / "downloads").glob("*/*.csv"))
    if len(found) != 1:
        raise ValueError(f"Expected one leaderboard CSV for {run}")
    return found[0]


inputs = {
    "before-d6.csv": single_csv("20260926T100000167581Z-92413ba1"),
    "after-d6.csv": single_csv("20260926T100306712074Z-6e0eb8fc"),
    "before-ma.csv": single_csv("20260926T095947592805Z-7d9bacdc"),
    "after-ma.csv": single_csv("20260926T100322667181Z-7a237ebe"),
    "before.csv": root
    / "artifacts/submissions/regression/20260924T054343781258Z-db4d3cf3/predictions.csv",
    "after.csv": root
    / "artifacts/submissions/regression/20260926T090548323924Z-46188564/predictions.csv",
    "receipt.json": root
    / "artifacts/submission-attempts/regression/20260926T100008522917Z-0ab7c954/result.json",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit reflected D6-only Live leaderboard change after robust-bound submission",
    "--config",
    str(root / "configs/regression/d6-robust-live-result-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_robust_live_result.py"],
    cwd=root,
    check=True,
)
