"""Replay the generic classification preflight against the failed D6 candidate."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
consensus = runs / "20260918T085712373599Z-8829d528"
support = runs / "20260918T085333214169Z-892612f4"
inputs = {
    "decision-test.parquet": consensus / "outputs/test.parquet",
    "decision-oof.parquet": consensus / "outputs/oof.parquet",
    "member-test.parquet": consensus / "inputs/split-bag-test.parquet",
    "member-oof.parquet": consensus / "inputs/split-bag-oof.parquet",
    "compounds.parquet": support / "inputs/compounds.parquet",
    "measurements.parquet": support / "inputs/measurements.parquet",
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
    "Replay the generic classification candidate preflight against failed D6 support6",
    "--config",
    str(root / "configs/analysis/d6-support6-preflight-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_candidate_preflight.py"],
    cwd=root,
    check=True,
)
