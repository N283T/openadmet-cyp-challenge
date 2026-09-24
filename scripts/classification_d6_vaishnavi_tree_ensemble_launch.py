"""Launch the public shallow-tree CYP2D6 reproduction."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "classification-labels.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "test-mapping.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "current-prediction-bank.parquet": root
    / (
        "artifacts/runs/classification/20260915T121009912691Z-73198c15/outputs/"
        "prediction-bank.parquet"
    ),
}
for path in inputs.values():
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
        "Reproduce the vaishnavi53 shallow RF/ExtraTrees/LightGBM CYP2D6 ensemble "
        "with and without the documented 13-binary alert approximation on UMAP "
        "and scaffold folds"
    ),
    "--config",
    str(root / "configs/classification/vaishnavi-tree-ensemble-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_vaishnavi_tree_ensemble.py"],
    cwd=root,
    check=True,
)
