"""Launch frozen CheMeleon embedding extraction for BioTransformer products."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root
    / (
        "artifacts/runs/classification/20260914T060605082765Z-9e3ff5cd/outputs/"
        "biotransformer-products.parquet"
    ),
    root / "artifacts/runs/shared/20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt",
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
    "Extract frozen CheMeleon embeddings for parents and one-step BioTransformer products",
    "--config",
    str(root / "configs/classification-round/biotransformer-product-embedding-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/biotransformer_product_embedding.py"],
    cwd=root,
    check=True,
)
