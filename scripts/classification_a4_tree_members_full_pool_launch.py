"""Record complete-pool inference for the submitted A4 tree members."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
lgbm = runs / "classification/20260912T115220209511Z-76facc71/outputs"
config = root / "configs/classification/a4-tree-members-full-pool-20260926.json"
inputs = {
    "compounds.parquet": runs / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "class-embeddings.parquet": runs
    / "shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet",
    "pxr-features.parquet": runs
    / "classification/20260913T152152397641Z-5177d5d9/outputs/pxr-classical-features.parquet",
    "admetica-features.parquet": runs
    / "regression/20260911T145924198895Z-373de784/outputs/score-features.parquet",
    "lgbm-oof.parquet": lgbm / "oof.parquet",
    "pxr-oof.parquet": runs / "classification/20260913T152746349812Z-c510a24e/outputs/oof.parquet",
}
for path in sorted(
    (runs / "shared/20260909T132918829974Z-bc81af58/outputs").glob("mordred-*.parquet")
):
    inputs[path.name] = path
for fold in range(5):
    inputs[f"external-class-mordred-fold{fold}.txt"] = (
        lgbm / f"external_class_mordred_balanced-CYP3A4_is_TDI-{fold}.txt"
    )
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
    "Extend exact A4 external-class/Mordred and PXR/Admetica trees to the complete UMAP pool",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
completed = subprocess.run(
    command + ["--", "python", "scripts/classification_a4_tree_members_full_pool.py"],
    cwd=root,
    check=False,
)
if completed.returncode:
    raise SystemExit(f"A4 tree-member run failed with exit code {completed.returncode}")
