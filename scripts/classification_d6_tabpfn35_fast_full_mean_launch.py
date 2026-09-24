"""Launch the matched TabPFN 3.5 fast/full mean audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
inputs = {
    "base-oof.parquet": classification / "20260918T063800192024Z-c1058115/outputs/oof.parquet",
    "deployed-test.parquet": classification
    / "20260915T125939829893Z-1d921ffc/outputs/test-scores.parquet",
    "fast-oof.parquet": classification / "20260918T063048677015Z-e0583fe4/outputs/oof.parquet",
    "fast-test.parquet": classification / "20260918T063048677015Z-e0583fe4/outputs/test.parquet",
    "full-oof.parquet": classification / "20260918T063852992056Z-387d4e0e/outputs/oof.parquet",
    "full-test.parquet": classification / "20260918T063852992056Z-387d4e0e/outputs/test.parquet",
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
    "Audit matched fast and full TabPFN 3.5 mean for CYP2D6 classification",
    "--config",
    str(root / "configs/classification/d6-tabpfn35-fast-full-mean-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_tabpfn35_fast_full_mean.py"],
    cwd=root,
    check=True,
)
