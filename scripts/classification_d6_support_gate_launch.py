"""Launch recorded support-aware CYP2D6 routing analysis."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    full = root / path
    assert full.is_file(), full
    paths.append(full)


shared = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
for filename in [
    "compounds.parquet",
    "measurements.parquet",
    "test_mapping.parquet",
]:
    add(f"{shared}/{filename}")
add("artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet")

diverse = "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs"
fixed = "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs"
fp = "artifacts/runs/classification/20260913T122320230702Z-d08c8ee5/outputs"
for path in [
    f"{diverse}/oof.parquet",
    f"{diverse}/test.parquet",
    f"{fixed}/oof.parquet",
    f"{fp}/oof.parquet",
    f"{fp}/test.parquet",
]:
    add(path)
for family in ["complement_nb_a1", "rbf_svm_balanced_c1"]:
    for fold in range(5):
        add(f"{fixed}/{family}-fold{fold}.joblib")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Test fold-excluded support-aware routing among saved D6 fingerprint, "
        "descriptor and Chemprop predictions without new model fitting"
    ),
    "--config",
    str(root / "configs/classification-round/d6-support-gate-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_support_gate.py"],
    cwd=root,
    check=True,
)
