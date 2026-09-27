"""Record reflected Live metrics for the combined calibration submission."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
api = root / "artifacts/api/regression"


def single_csv(run: str) -> Path:
    found = list((api / run / "downloads").glob("*/*.csv"))
    if len(found) != 1:
        raise ValueError(f"Expected one leaderboard CSV for {run}")
    return found[0]


runs = {
    "MA": ("20260927T051400483984Z-c461cf12", "20260927T052855223772Z-2c84697f"),
    "1A2": ("20260927T051406426900Z-be611869", "20260927T052855224585Z-a1d89027"),
    "2C9": ("20260927T051410837796Z-ff4b6690", "20260927T052855223778Z-84cd1ae5"),
    "2D6": ("20260927T051414821496Z-2b02f97b", "20260927T052855223778Z-38c108eb"),
    "3A4": ("20260927T051419350141Z-6d6338f6", "20260927T052855223785Z-f254c1c5"),
}
inputs = {
    "before.csv": root
    / "artifacts/submissions/regression/20260926T090548323924Z-46188564/predictions.csv",
    "after.csv": root
    / "artifacts/submissions/regression/20260926T140601457029Z-bd3a1e40/predictions.csv",
    "receipt.json": root
    / "artifacts/submission-attempts/regression/20260927T051423190683Z-117605b9/result.json",
}
for board, (old_run, new_run) in runs.items():
    inputs[f"before-{board}.csv"] = single_csv(old_run)
    inputs[f"after-{board}.csv"] = single_csv(new_run)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit reflected Live endpoint and macro results of combined affine calibration",
    "--config",
    str(root / "configs/regression/calibration-live-result-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_calibration_live_result.py"],
    cwd=root,
    check=True,
)
