"""Record a fixed weak-probability adjustment against matched regression OOF."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "classifier-oof.parquet": (
        runs / "regression/20260926T014119877434Z-bf53d48d/outputs/oof.parquet"
    ),
    "parent-base-oof.parquet": (
        runs / "regression/20260918T012518053809Z-5c46e72d/outputs/jeremy_anchor_all-oof.parquet"
    ),
    "scaffold-base-oof.parquet": (
        runs / "regression/20260918T021832961989Z-c89b25fa/outputs/jeremy_anchor_all-oof.parquet"
    ),
    "parent-assignments.parquet": (
        runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "scaffold-assignments.parquet": (
        runs / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet"
    ),
    "measurements.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Audit fixed weak-probability corrections on parent and scaffold OOF",
    "--config", str(root / "configs/regression/d6-weak-probability-correction-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_weak_probability_correction.py"],
    cwd=root,
    check=True,
)
