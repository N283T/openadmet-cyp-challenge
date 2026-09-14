"""Launch the recorded A4 FDA source-support gate experiment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / (
        "artifacts/runs/classification/20260913T143057530832Z-47a95438/outputs/"
        "fda-cyp3a4-tdi-records.parquet"
    ),
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/classification/20260913T143341518156Z-5463a4dc/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T143341518156Z-5463a4dc/outputs/test.parquet",
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
        "Test whether Morgan similarity support gating makes the independent FDA "
        "CYP3A4 TDI LR block directionally stable across fixed UMAP and scaffold "
        "validation"
    ),
    "--config",
    str(root / "configs/classification-round/a4-fda-support-gate-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_fda_support_gate.py"],
    cwd=root,
    check=True,
)
