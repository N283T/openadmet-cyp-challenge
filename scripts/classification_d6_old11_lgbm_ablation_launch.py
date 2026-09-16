"""Launch exact fixed-policy ablation of D6 old11 LightGBM members."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
shared = root / "artifacts/runs/shared"
historical = classification / "20260912T121146867140Z-b0a05986"
fixed = classification / "20260913T123045515718Z-0df22d54/outputs"
inputs = [
    historical / "outputs/ensemble-oof.parquet",
    historical / "outputs/ensemble-test.parquet",
    historical / "inputs/000-oof.parquet",
    historical / "inputs/002-fold-test-probabilities.parquet",
    historical / "inputs/003-oof-probabilities.parquet",
    historical / "inputs/005-test-fold-probabilities.parquet",
    historical / "inputs/006-classification-oof.parquet",
    historical / "inputs/008-fold-test-probabilities.parquet",
    fixed / "oof.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    *[fixed / f"logistic_c0.05-fold{fold}.joblib" for fold in range(5)],
    *[fixed / f"complement_nb_a1-fold{fold}.joblib" for fold in range(5)],
    *[fixed / f"rbf_svm_c1-fold{fold}.joblib" for fold in range(5)],
]
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
    "Ablate basic Morgan and Mordred LightGBM members from the fixed D6 legacy hierarchy",
    "--config",
    str(root / "configs/classification-round/d6-old11-lgbm-ablation-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_old11_lgbm_ablation.py"],
    cwd=root,
    check=True,
)
