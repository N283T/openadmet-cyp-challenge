"""Record the September 15 control composition on complete UMAP folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-sep15-control-full-pool-v1.json"
inputs = {
    "lgbm-all-pool.parquet": runs / "20260924T035725764059Z-9b7d6d23/outputs/all-pool-oof.parquet",
    "pfn-pca-all-pool.parquet": runs
    / "20260924T034908336839Z-ec9287a1/outputs/all-pool-oof.parquet",
    "pfn-raw-all-pool.parquet": runs
    / "20260924T035915695580Z-7dd31f09/outputs/all-pool-oof.parquet",
    "graph-all-pool.parquet": runs / "20260924T040130129624Z-558f74b9/outputs/all-pool-oof.parquet",
    "fixed-all-pool.parquet": runs / "20260924T040107573312Z-887ed73b/outputs/all-pool-oof.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {declared[key] for key in ["lgbm", "pfn_pca", "pfn_raw", "graph", "fixed"]}:
    raise ValueError("Config and input aliases differ")
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
    "Recompose September 15 D6 control on complete UMAP held-fold pools",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_sep15_control_full_pool.py"],
    cwd=root,
    check=True,
)
