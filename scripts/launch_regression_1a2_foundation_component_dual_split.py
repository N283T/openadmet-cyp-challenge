"""Audit selected CYP1A2 atomic additions on both existing split maps."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
primary_lgbm = runs / "20260918T012518053809Z-5c46e72d/outputs"
scaffold_lgbm = runs / "20260918T021832961989Z-c89b25fa/outputs"
primary_pfn = runs / "20260918T030037649196Z-3101f488/outputs"
scaffold_pfn = runs / "20260918T033058771796Z-79d216a4/outputs"
inputs = {
    "parent-all_pool_lgbm.parquet": primary_lgbm / "all_pool_top512-oof.parquet",
    "parent-dense_pool_lgbm.parquet": primary_lgbm / "dense_pool_top256-oof.parquet",
    "parent-all_pool_tabpfn.parquet": (
        primary_pfn / "all_pool_top512__tabpfn35fast-CYP1A2-oof.parquet"
    ),
    "parent-jeremy_anchor_lgbm.parquet": primary_lgbm / "jeremy_anchor_all-oof.parquet",
    "scaffold-all_pool_lgbm.parquet": scaffold_lgbm / "all_pool_top512-oof.parquet",
    "scaffold-dense_pool_lgbm.parquet": scaffold_lgbm / "dense_pool_top256-oof.parquet",
    "scaffold-all_pool_tabpfn.parquet": (
        scaffold_pfn / "all_pool_top512__tabpfn35fast-CYP1A2-oof.parquet"
    ),
    "scaffold-jeremy_anchor_lgbm.parquet": scaffold_lgbm / "jeremy_anchor_all-oof.parquet",
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Compare CYP1A2 foundation and LightGBM atomic members on both splits",
    "--config", str(root / "configs/regression/1a2-foundation-component-dual-split-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_1a2_foundation_component_dual_split.py"],
    cwd=root,
    check=True,
)
