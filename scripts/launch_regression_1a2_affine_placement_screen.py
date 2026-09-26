"""Record affine sensitivity on the current 1A2 parent bank and dual-split proxy."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg_path = root / "configs/regression/1a2-affine-placement-screen-v1.json"
cfg = json.loads(cfg_path.read_text())
inputs = {
    "parent-oof.parquet": root
    / "artifacts/runs/regression"
    / cfg["parent_oof_run"]
    / "outputs/1A2/add_dense_and_fast_pool/oof.parquet",
    "proxy-oof.parquet": root
    / "artifacts/runs/regression"
    / cfg["proxy_oof_run"]
    / "outputs/oof.parquet",
    "current.csv": root
    / "artifacts/submissions/regression"
    / cfg["current_submission_id"]
    / "predictions.csv",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Screen 1A2 location and spread on current parent bank and dual-split three-member proxy",
    "--config",
    str(cfg_path),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_1a2_affine_placement_screen.py"],
    cwd=root,
    check=True,
)
