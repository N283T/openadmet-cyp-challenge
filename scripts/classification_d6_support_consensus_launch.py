"""Launch equal-vote support-relevant D6 consensus audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = {
    "split-bag-oof.parquet": runs / "20260918T085533538849Z-fb827f03/outputs/oof.parquet",
    "split-bag-test.parquet": runs / "20260918T085533538849Z-fb827f03/outputs/test.parquet",
    "random-oof-support.parquet": runs
    / "20260918T085333214169Z-892612f4/outputs/random-oof-support.parquet",
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
    "Audit equal top-call consensus across complete CYP2D6 TabPFN split grids",
    "--config",
    str(root / "configs/classification/d6-support-consensus-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_support_consensus.py"],
    cwd=root,
    check=True,
)
