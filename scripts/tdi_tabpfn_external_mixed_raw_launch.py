"""Launch raw-600 external-mixed TabPFN on D6."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
prep = runs / "classification/20260912T115018743045Z-269e598d/outputs"
inputs = [
    prep / "classification-labels.parquet",
    prep / "umap-assignments.parquet",
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    runs
    / "shared/20260910T121624900089Z-271684c2/inputs"
    / "004-tabpfn-v2.6-classifier-v2.6_default.ckpt",
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
    "Compare raw scaled 600-dimensional external-mixed embeddings with PCA256 for D6 TabPFN",
    "--config",
    str(root / "configs/classification-tabpfn/external-mixed-raw600-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/tdi_tabpfn_frozen.py"], cwd=root, check=True)
