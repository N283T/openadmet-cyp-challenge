"""Evaluate target-specific v7 swaps and foundation-readout additions."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs_dir = runs / "20260918T015803413341Z-7286be9d/inputs"
lgbm = runs / "20260918T012518053809Z-5c46e72d/outputs"
foundation = runs / "20260918T030037649196Z-3101f488/outputs"
config = root / "configs/regression/flat-bank-foundation-swap-v8.json"
inputs = {path.name: path for path in inputs_dir.iterdir() if path.is_file()}
inputs.update(
    {
        "lgbm-jeremy-all-log2fc.parquet": lgbm / "jeremy_all_log2fc-oof.parquet",
        "lgbm-jeremy-anchor-all-log2fc.parquet": lgbm / "jeremy_anchor_all_log2fc-oof.parquet",
        "lgbm-dense-pool-top512.parquet": lgbm / "dense_pool_top512-oof.parquet",
        "foundation-fast-all512-1a2.parquet": foundation
        / "all_pool_top512__tabpfn35fast-CYP1A2-oof.parquet",
        "foundation-fast-jeremy-anchor-1a2.parquet": foundation
        / "jeremy_anchor_all__tabpfn35fast-CYP1A2-oof.parquet",
        "foundation-fast-dense512-2c9.parquet": foundation
        / "dense_pool_top512__tabpfn35fast-CYP2C9-oof.parquet",
        "foundation-icl-dense512-2c9.parquet": foundation
        / "dense_pool_top512__tabicl-CYP2C9-oof.parquet",
        "foundation-fast-dense512-3a4.parquet": foundation
        / "dense_pool_top512__tabpfn35fast-CYP3A4-oof.parquet",
        "foundation-icl-dense512-3a4.parquet": foundation
        / "dense_pool_top512__tabicl-CYP3A4-oof.parquet",
    }
)
declared = json.loads(config.read_text())
needed = {declared["scorer"], declared["current_control"]}
for pool in declared["pools"].values():
    needed.update(file for _, file in pool["components"])
missing = needed - set(inputs)
if missing:
    raise FileNotFoundError(sorted(missing))

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare reviewed v7 banks with LGBM swaps and foundation-readout additions",
    "--config",
    str(config),
]
for alias in sorted(needed):
    command.extend(["--input", f"{alias}={inputs[alias]}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_flat_bank_review.py"],
    cwd=root,
    check=True,
)
