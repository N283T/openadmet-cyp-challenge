"""Launch fold-safe activity scoring of BioTransformer products."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
products = root / "artifacts/runs/classification/20260914T060605082765Z-9e3ff5cd/outputs"
inputs = [
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    shared / "20260910T100337258849Z-927eac1c/outputs/targets.parquet",
    shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    products / "biotransformer-products.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
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
        "Score one-step BioTransformer products using fold-safe D6/A4 direct-pIC50 "
        "and primary-log2fc Morgan regressors"
    ),
    "--config",
    str(root / "configs/classification-round/biotransformer-product-activity-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/biotransformer_product_activity.py"],
    cwd=root,
    check=True,
)
