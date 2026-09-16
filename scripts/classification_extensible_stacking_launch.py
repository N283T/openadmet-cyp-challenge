"""Launch the prespecified extensible classification stacking comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
historical = classification / "20260912T121146867140Z-b0a05986"
fixed = classification / "20260913T123045515718Z-0df22d54/outputs"
inputs = [
    historical / "outputs/ensemble-oof.parquet",
    historical / "inputs/000-oof.parquet",
    historical / "inputs/003-oof-probabilities.parquet",
    historical / "inputs/006-classification-oof.parquet",
    fixed / "oof.parquet",
    classification / "20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    classification / "20260913T152746349812Z-c510a24e/outputs/oof.parquet",
    classification / "20260914T063614433212Z-c0dc3f32/outputs/oof.parquet",
    classification / "20260915T065610180920Z-f532943e/outputs/oof.parquet",
    classification / "20260914T051014577658Z-7b1b64c2/outputs/oof.parquet",
    classification / "20260915T073340818290Z-d781ddd2/outputs/oof.parquet",
    classification / "20260915T112944344816Z-69abecf5/outputs/oof-probabilities.parquet",
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
    "Compare extensible ridge, convex-simplex, pairwise, and expert-set stackers",
    "--config",
    str(root / "configs/classification-round/extensible-stacking-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_extensible_stacking.py"],
    cwd=root,
    check=True,
)
