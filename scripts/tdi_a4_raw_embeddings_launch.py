"""Launch matched raw600 A4 TabICL and TabPFN controls."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
old = runs / "classification/20260912T115242271820Z-51023dac/inputs"
inputs = [
    old / "000-classification-labels.parquet",
    old / "001-umap-assignments.parquet",
    old / "002-test_mapping.parquet",
    old / "003-004-tabpfn-v2.6-classifier-v2.6_default.ckpt",
    runs / "shared/20260912T134236392057Z-8c08c1de/outputs/tabicl-classifier-v2-20260212.ckpt",
    runs / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    runs / "shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet",
    runs / "shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet",
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
    "Fit matched raw600 A4 TabICL8 and TabPFN8 outer and inner controls",
    "--config",
    str(root / "configs/classification-tabicl/a4-raw-embeddings-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/tdi_jeremy_readout.py"], cwd=root, check=True)
