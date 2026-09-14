"""Launch diverse non-LightGBM CYP2D6 rankers with complete input lineage."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    full = root / path
    assert full.is_file(), full
    paths.append(full)


shared = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
for filename in ["compounds.parquet", "measurements.parquet", "test_mapping.parquet"]:
    add(f"{shared}/{filename}")
add("artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet")
add("artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet")
add("artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet")
for path in sorted(
    (root / "artifacts/runs/shared/20260909T132918829974Z-bc81af58/outputs").glob(
        "mordred-*.parquet"
    )
):
    add(path.relative_to(root))
add("artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet")
add("artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/ensemble-test.parquet")
add("artifacts/runs/classification/20260913T101151963949Z-892b1e27/outputs/summary.json")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Compare diverse CatBoost, ExtraTrees, logistic and Tanimoto CYP2D6 "
        "rankers on fixed UMAP folds"
    ),
    "--config",
    str(root / "configs/classification-round/d6-diverse-models-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_diverse_models.py"],
    cwd=root,
    check=True,
)
