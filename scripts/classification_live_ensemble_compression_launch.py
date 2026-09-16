"""Launch expanded-member compression of the current classification model."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
b0 = classification / "20260912T121146867140Z-b0a05986"
fixed = classification / "20260913T123045515718Z-0df22d54"
inputs = [
    classification / "20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    classification / "20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
    classification / "20260913T152746349812Z-c510a24e/outputs/oof.parquet",
    classification / "20260913T152746349812Z-c510a24e/outputs/test.parquet",
    b0 / "outputs/ensemble-oof.parquet",
    b0 / "outputs/ensemble-test.parquet",
    b0 / "inputs/000-oof.parquet",
    b0 / "inputs/002-fold-test-probabilities.parquet",
    b0 / "inputs/003-oof-probabilities.parquet",
    b0 / "inputs/005-test-fold-probabilities.parquet",
    b0 / "inputs/006-classification-oof.parquet",
    b0 / "inputs/008-fold-test-probabilities.parquet",
    fixed / "outputs/oof.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
]
inputs.extend(fixed / "outputs" / f"logistic_c0.05-fold{fold}.joblib" for fold in range(5))
inputs.extend(fixed / "outputs" / f"complement_nb_a1-fold{fold}.joblib" for fold in range(5))
for path in inputs:
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
    "Distill current D6 and A4 classifiers over their expanded member sets",
    "--config",
    str(root / "configs/classification-round/live-ensemble-compression-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_live_ensemble_compression.py"],
    cwd=root,
    check=True,
)
