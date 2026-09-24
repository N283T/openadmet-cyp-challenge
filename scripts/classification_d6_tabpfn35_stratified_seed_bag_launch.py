"""Launch the support-selected five-seed D6 TabPFN bag audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = {
    "primary-oof.parquet": runs / "20260918T092956435999Z-277c5b65/outputs/oof.parquet",
    "primary-test.parquet": runs / "20260918T092956435999Z-277c5b65/outputs/test.parquet",
    "followup-oof.parquet": runs / "20260918T093503503912Z-03186dea/outputs/oof.parquet",
    "followup-test.parquet": runs / "20260918T093503503912Z-03186dea/outputs/test.parquet",
    "deployed-oof.parquet": runs / "20260918T063800192024Z-c1058115/outputs/oof.parquet",
    "deployed-test.parquet": runs
    / "20260915T125939829893Z-1d921ffc/outputs/test-scores.parquet",
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
    "Audit equal five-seed support-selected stratified-random CYP2D6 TabPFN bags",
    "--config",
    str(root / "configs/classification/d6-tabpfn35-stratified-seed-bag-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_tabpfn35_stratified_seed_bag.py"],
    cwd=root,
    check=True,
)
