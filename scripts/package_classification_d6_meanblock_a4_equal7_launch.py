"""Launch the recorded packaging run for the two selected classification candidates."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
historical = classification / "20260912T121146867140Z-b0a05986/inputs"
fixed = classification / "20260913T123045515718Z-0df22d54/outputs"

inputs: list[Path] = []


def add(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    name = f"{len(inputs):03d}-{path.name}"
    inputs.append(path)
    return name


cfg = {
    "lgbm_test": add(historical / "002-fold-test-probabilities.parquet"),
    "pfn_test": add(historical / "005-test-fold-probabilities.parquet"),
    "graph_test": add(historical / "008-fold-test-probabilities.parquet"),
    "raw600_test": add(
        classification
        / "20260915T112944344816Z-69abecf5/outputs/test-probabilities.parquet"
    ),
    "compounds": add(shared / "compounds.parquet"),
    "test_mapping": add(shared / "test_mapping.parquet"),
    "aliases": add(shared / "aliases.parquet"),
    "a4_test": add(classification / "20260915T125001009276Z-d5adff22/outputs/test.parquet"),
    "reference": add(root / "data/raw/cyp-challenge-TEST-BLINDED.csv"),
    "logistic_models": [add(fixed / f"logistic_c0.05-fold{fold}.joblib") for fold in range(5)],
    "rbf_models": [add(fixed / f"rbf_svm_c1-fold{fold}.joblib") for fold in range(5)],
    "complement_models": [
        add(fixed / f"complement_nb_a1-fold{fold}.joblib") for fold in range(5)
    ],
    "d6_test_calls": 49,
    "a4_model": "exhaustive_mcc::components_7",
    "a4_test_calls": 249,
    "selection": (
        "Package two prespecified candidates: the D6 PCA256/raw600 equal embedding block "
        "inside the full-data ten-component selection, and the A4 fixed equal-seven subset."
    ),
}
config = root / "configs/submission-classification/d6-meanblock-a4-equal7-v1.json"
if not config.is_file():
    raise FileNotFoundError(config)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Freeze D6 mean-block and A4 equal-seven classification submission",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/package_classification_d6_meanblock_a4_equal7.py"],
    cwd=root,
    check=True,
)
