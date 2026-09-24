"""Record the frozen raw600 TabPFN control atom on complete UMAP folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260915T112944344816Z-69abecf5"
config = root / "configs/analysis/d6-control-pfn-raw600-all-pool-v1.json"
cfg = json.loads(config.read_text())
inputs = {
    "labels.parquet": original / "inputs/000-classification-labels.parquet",
    "assignments.parquet": original / "inputs/001-umap-assignments.parquet",
    "external-embeddings.parquet": original / "inputs/004-external-embeddings.parquet",
    "tabpfn-v2.6.ckpt": original / "inputs/003-004-tabpfn-v2.6-classifier-v2.6_default.ckpt",
}
for fold, alias in enumerate(cfg["preprocess_aliases"]):
    name = f"external_mixed_raw600-CYP2D6_is_TDI-fold{fold}-pca.joblib"
    inputs[alias] = original / "outputs" / name
for path in [config, *inputs.values()]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi", "run", "cyp", "run", "--track", "classification", "--purpose",
    "Refit September 15 raw600 TabPFN control atom on complete UMAP folds",
    "--config", str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_legacy_pfn_all_pool.py"],
    cwd=root,
    check=True,
)
