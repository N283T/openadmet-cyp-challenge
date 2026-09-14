"""Launch the recorded deployed A4 Caruana member-correlation audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification/20260912T141520708830Z-97df4c71"
inputs = [
    run / "outputs/selections.json",
    run / "inputs/000-000-000-oof.parquet",
    run / "inputs/002-002-002-fold-test-probabilities.parquet",
    run / "inputs/003-003-003-oof-probabilities.parquet",
    run / "inputs/005-005-005-test-fold-probabilities.parquet",
    run / "inputs/009-009-outer-oof.parquet",
    run / "inputs/011-011-fold-test-probabilities.parquet",
    run / "inputs/015-outer-oof.parquet",
    run / "inputs/017-fold-test-probabilities.parquet",
]
for path in inputs:
    assert path.is_file(), path
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Measure OOF and blind-test prediction dependence, call overlap, and standalone "
        "A4 performance for the six deployed Caruana-AUROC members"
    ),
    "--config",
    str(root / "configs/classification-round/caruana-member-correlation-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_caruana_member_correlation.py"],
    cwd=root,
    check=True,
)
