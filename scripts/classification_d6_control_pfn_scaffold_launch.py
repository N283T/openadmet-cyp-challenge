"""Record scaffold refits of the September 15 PCA and raw TabPFN atoms."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("view", choices=["pca", "raw"])
view = parser.parse_args().view
root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260912T115242271820Z-51023dac"
config = root / f"configs/analysis/d6-control-pfn-scaffold-{view}-v1.json"
inputs = {
    "labels.parquet": original / "inputs/000-classification-labels.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "external-embeddings.parquet": original / "inputs/005-external-embeddings.parquet",
    "tabpfn-v2.6.ckpt": original / "inputs/003-004-tabpfn-v2.6-classifier-v2.6_default.ckpt",
}
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
    f"Refit September 15 {view} TabPFN control atom on complete scaffold folds",
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
