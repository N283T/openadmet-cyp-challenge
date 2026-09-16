"""Launch current-live A4 blend audit for CheMeleon parent/product embeddings."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = [
    runs / "20260913T153156161201Z-1954e23e/inputs/000-member-oof.parquet",
    runs / "20260913T153156161201Z-1954e23e/inputs/002-oof.parquet",
    runs / "20260913T153156161201Z-1954e23e/outputs/test-blends.parquet",
    runs / "20260915T070440285192Z-321d27f7/outputs/oof.parquet",
    runs / "20260915T070440285192Z-321d27f7/outputs/test.parquet",
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
    "Audit fixed CheMeleon parent/product embedding additions to current A4 rank",
    "--config",
    str(
        root
        / "configs/classification-round/a4-biotransformer-embedding-live-blend-v1.json"
    ),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command
    + ["--", "python", "scripts/classification_a4_biotransformer_embedding_blend.py"],
    cwd=root,
    check=True,
)
