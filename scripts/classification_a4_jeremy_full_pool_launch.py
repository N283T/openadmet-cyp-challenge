"""Record all-molecule inference from the saved Jeremy600 A4 fold models."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
pca = runs / "classification/20260912T140401138774Z-0e4bcc8d/outputs"
raw = runs / "classification/20260915T123533598879Z-a459b098/outputs"
weights = runs / "classification/20260912T140401138774Z-0e4bcc8d/inputs"
config = root / "configs/classification/a4-jeremy-full-pool-20260926.json"
cfg = json.loads(config.read_text())
inputs = {
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "jeremy-features.parquet": runs
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "pca-oof.parquet": pca / "outer-oof.parquet",
    "raw-oof.parquet": raw / "outer-oof.parquet",
    "pfn-weights.ckpt": weights / "003-003-004-tabpfn-v2.6-classifier-v2.6_default.ckpt",
    "icl-weights.ckpt": weights / "004-tabicl-classifier-v2-20260212.ckpt",
}
for view in cfg["views"]:
    phase, backend, method = view["phase"], view["backend"], view["method"]
    source = pca if phase == "pca" else raw
    for fold in range(5):
        original = f"{method}-CYP3A4_is_TDI-fold{fold}-outer"
        alias = f"{phase}-{backend}-fold{fold}"
        inputs[f"{alias}-prep.joblib"] = source / f"{original}-pca.joblib"
        extension = "tabicl_fit" if backend == "tabicl" else "tabpfn_fit"
        inputs[f"{alias}.fit"] = source / f"{original}.{extension}"
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
    "Extend all four saved Jeremy A4 readouts to every UMAP held-fold molecule",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
completed = subprocess.run(
    command + ["--", "python", "scripts/classification_a4_jeremy_full_pool.py"],
    cwd=root,
    check=False,
)
if completed.returncode:
    raise SystemExit(f"Jeremy A4 full-pool run failed with exit code {completed.returncode}")
