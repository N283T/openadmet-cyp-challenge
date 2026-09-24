"""Record the frozen legacy TabPFN AP atom on the complete UMAP pool."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260912T115242271820Z-51023dac"
config = root / "configs/analysis/d6-ap-legacy-pfn-all-pool-v1.json"
cfg = json.loads(config.read_text())
inputs = {
    "labels.parquet": original / "inputs/000-classification-labels.parquet",
    "assignments.parquet": original / "inputs/001-umap-assignments.parquet",
    "external-embeddings.parquet": original / "inputs/005-external-embeddings.parquet",
    "tabpfn-v2.6.ckpt": original / "inputs/003-004-tabpfn-v2.6-classifier-v2.6_default.ckpt",
}
for fold in range(5):
    prefix = f"external_mixed_pca256-CYP2D6_is_TDI-fold{fold}"
    inputs[cfg["preprocess_aliases"][fold]] = original / "outputs" / f"{prefix}-pca.joblib"
for path in [config, *inputs.values()]:
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
    "Refit exact legacy TabPFN 2.6 AP recipe on complete held-fold pool",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_legacy_pfn_all_pool.py"],
    cwd=root,
    check=True,
)
