"""Launch cross-fit peer-pool and one-for-one replacement audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
shared = root / "artifacts/runs/shared"
fixed = classification / "20260913T123045515718Z-0df22d54/outputs"
historical = classification / "20260912T121146867140Z-b0a05986/inputs"
inputs = [
    classification / "20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    classification / "20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
    classification / "20260913T152746349812Z-c510a24e/outputs/oof.parquet",
    classification / "20260913T152746349812Z-c510a24e/outputs/test.parquet",
    classification / "20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet",
    classification / "20260912T121146867140Z-b0a05986/outputs/ensemble-test.parquet",
    fixed / "oof.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    *[fixed / f"logistic_c0.05-fold{fold}.joblib" for fold in range(5)],
    *[fixed / f"complement_nb_a1-fold{fold}.joblib" for fold in range(5)],
    classification / "20260914T063614433212Z-c0dc3f32/outputs/oof.parquet",
    classification / "20260914T063614433212Z-c0dc3f32/outputs/test.parquet",
    classification / "20260915T065610180920Z-f532943e/outputs/oof.parquet",
    classification / "20260915T065610180920Z-f532943e/outputs/test.parquet",
    classification / "20260914T051014577658Z-7b1b64c2/outputs/oof.parquet",
    classification / "20260914T051014577658Z-7b1b64c2/outputs/test.parquet",
    classification / "20260915T073340818290Z-d781ddd2/outputs/oof.parquet",
    classification / "20260915T073340818290Z-d781ddd2/outputs/test.parquet",
    *[fixed / f"rbf_svm_c1-fold{fold}.joblib" for fold in range(5)],
    historical / "000-oof.parquet",
    historical / "002-fold-test-probabilities.parquet",
    historical / "003-oof-probabilities.parquet",
    historical / "005-test-fold-probabilities.parquet",
    historical / "006-classification-oof.parquet",
    historical / "008-fold-test-probabilities.parquet",
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
    (
        "Cross-fit previously rejected representation families as Caruana peers "
        "and explicit replacements"
    ),
    "--config",
    str(root / "configs/classification-round/live-peer-replacement-crossfit-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_live_peer_replacement_crossfit.py"],
    cwd=root,
    check=True,
)
