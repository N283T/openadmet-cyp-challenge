"""Launch TDI readouts over frozen parent/product CheMeleon embeddings."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
embedding = root / "artifacts/runs/classification/20260915T070040196804Z-517e5bfc/outputs"
inputs = [
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    embedding / "product-embeddings.parquet",
    embedding / "product-embedding-scalars.parquet",
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
    "Evaluate frozen CheMeleon parent and BioTransformer-product embedding differences",
    "--config",
    str(
        root
        / "configs/classification-round/biotransformer-product-embedding-readout-v1.json"
    ),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_biotransformer_embedding_readout.py"],
    cwd=root,
    check=True,
)
