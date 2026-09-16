"""Launch the generic readout for one frozen representation family."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument(
    "family",
    choices=[
        "attnsom",
        "mace_off",
        "mechanism_combined",
        "mechanism_plus_cypmap",
        "biotransformer",
        "cypmap",
    ],
)
args = parser.parse_args()

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
assignment = (
    root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
)
if args.family == "attnsom":
    features = root / "artifacts/runs/classification/20260914T045809423993Z-c5c2c254/outputs"
    family_inputs = [features / "attnsom-features.parquet", features / "attnsom-embeddings.parquet"]
elif args.family == "mace_off":
    features = root / "artifacts/runs/classification/20260914T050016310677Z-e80256e2/outputs"
    family_inputs = [
        features / "mace-off-properties.parquet",
        features / "mace-off-embeddings.parquet",
    ]
elif args.family == "mechanism_combined":
    smartcyp = root / "artifacts/runs/classification/20260914T044537363003Z-5e7f9e02/outputs"
    attnsom = root / "artifacts/runs/classification/20260914T045809423993Z-c5c2c254/outputs"
    family_inputs = [
        smartcyp / "smartcyp-features.parquet",
        attnsom / "attnsom-features.parquet",
        attnsom / "attnsom-embeddings.parquet",
    ]
elif args.family == "mechanism_plus_cypmap":
    smartcyp = root / "artifacts/runs/classification/20260914T044537363003Z-5e7f9e02/outputs"
    attnsom = root / "artifacts/runs/classification/20260914T045809423993Z-c5c2c254/outputs"
    cypmap = root / "artifacts/runs/classification/20260914T062717289785Z-17d2a645/outputs"
    family_inputs = [
        smartcyp / "smartcyp-features.parquet",
        attnsom / "attnsom-features.parquet",
        attnsom / "attnsom-embeddings.parquet",
        cypmap / "cypmap-features.parquet",
        cypmap / "cypmap-embeddings.parquet",
    ]
elif args.family == "biotransformer":
    features = root / "artifacts/runs/classification/20260914T060605082765Z-9e3ff5cd/outputs"
    family_inputs = [
        features / "biotransformer-features.parquet",
        features / "biotransformer-product-fingerprints.parquet",
    ]
else:
    features = root / "artifacts/runs/classification/20260914T062717289785Z-17d2a645/outputs"
    family_inputs = [
        features / "cypmap-features.parquet",
        features / "cypmap-embeddings.parquet",
    ]

inputs = [
    shared / "compounds.parquet",
    shared / "measurements.parquet",
    shared / "test_mapping.parquet",
    assignment,
    *family_inputs,
]
config = root / f"configs/classification-round/{args.family.replace('_', '-')}-readout-v1.json"
for path in [config, *inputs]:
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
    f"Evaluate frozen {args.family} representations on fixed UMAP and scaffold TDI folds",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_embedding_readout.py"], cwd=root, check=True
)
