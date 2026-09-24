"""Record all seven submitted AP atoms on complete scaffold folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-ap-full7-scaffold-full-pool-v1.json"
inputs = {
    "extra-all-pool.parquet": runs / "20260924T033028487589Z-c42e6f27/outputs/all-pool-oof.parquet",
    "lgbm-all-pool.parquet": runs / "20260924T033142058002Z-48b841b3/outputs/all-pool-oof.parquet",
    "tabpfn-all-pool.parquet": runs
    / "20260924T033614762367Z-eaac6dd8/outputs/all-pool-oof.parquet",
    "mechanism-all-pool.parquet": runs
    / "20260924T034339472109Z-62dde800/outputs/all-pool-oof.parquet",
    "activity-all-pool.parquet": runs
    / "20260924T034433332695Z-95f78765/outputs/all-pool-oof.parquet",
    "legacy-pfn-all-pool.parquet": runs
    / "20260924T041032506702Z-b92b2044/outputs/all-pool-oof.parquet",
    "ap-summary.json": runs / "20260923T173655387004Z-97d71dce/outputs/summary.json",
    "ap-outer-details.parquet": runs
    / "20260923T171350543926Z-9081bf31/outputs/outer-details.parquet",
}
declared = json.loads(config.read_text())
needed = {declared["ap_summary"], declared["outer_details"]} | {
    declared[member["source"]] for member in declared["members"]
}
if set(inputs) != needed:
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
    "Compose all seven AP atoms refit on complete scaffold held-fold pools",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_top3_full_pool.py"],
    cwd=root,
    check=True,
)
