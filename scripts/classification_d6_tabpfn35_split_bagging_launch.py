"""Launch the D6 TabPFN 3.5 equal split-bagging audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = {
    "dual-oof.parquet": runs / "20260918T063048677015Z-e0583fe4/outputs/oof.parquet",
    "dual-test.parquet": runs / "20260918T063048677015Z-e0583fe4/outputs/test.parquet",
    "random-oof.parquet": runs / "20260918T084809611158Z-8b586046/outputs/oof.parquet",
    "random-test.parquet": runs / "20260918T084809611158Z-8b586046/outputs/test.parquet",
    "random-oof-support.parquet": runs
    / "20260918T085333214169Z-892612f4/outputs/random-oof-support.parquet",
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
    "Audit equal UMAP scaffold random split-bagging for standalone CYP2D6 TabPFN 3.5",
    "--config",
    str(root / "configs/classification/d6-tabpfn35-split-bagging-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_tabpfn35_split_bagging.py"],
    cwd=root,
    check=True,
)
