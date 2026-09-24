"""Launch the A4 equal-seven plus full TabPFN 3.5 audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
inputs = {
    "equal7-oof.parquet": classification
    / "20260915T125001009276Z-d5adff22/outputs/oof.parquet",
    "equal7-test.parquet": classification
    / "20260915T125001009276Z-d5adff22/outputs/test.parquet",
    "tabpfn35-oof.parquet": classification
    / "20260918T071212026715Z-5c3e72b5/outputs/oof.parquet",
    "tabpfn35-test.parquet": classification
    / "20260918T071212026715Z-5c3e72b5/outputs/test.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
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
    "Audit full TabPFN 3.5 all-pool additions to the scored A4 equal-seven rank",
    "--config",
    str(root / "configs/classification/a4-tabpfn35-full-equal7-blend-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_tabpfn35_equal7_blend.py"],
    cwd=root,
    check=True,
)
