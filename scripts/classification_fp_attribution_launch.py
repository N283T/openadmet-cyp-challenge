"""Launch recorded Morgan-bit attribution and nested ablation analysis."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    full = root / path
    assert full.is_file(), full
    paths.append(full)


shared = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
for filename in ["compounds.parquet", "measurements.parquet"]:
    add(f"{shared}/{filename}")
add("artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet")

models = "artifacts/runs/classification/20260912T115220209511Z-76facc71/outputs"
for endpoint in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
    for fold in range(5):
        add(f"{models}/morgan_plain-{endpoint}-{fold}.txt")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Attribute stable Morgan bits in the retained TDI model components and "
        "evaluate radius and top-bit ablations on fixed UMAP folds"
    ),
    "--config",
    str(root / "configs/classification-round/fp-attribution-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_fp_attribution.py"],
    cwd=root,
    check=True,
)
