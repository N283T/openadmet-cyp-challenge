"""Record exact historical A4 stack inference on all held-fold molecules."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = root / "artifacts/runs/shared/20260910T132141471698Z-343eb5ec"
baseline = json.loads((source / "config.json").read_text())
config = root / "configs/classification/a4-historical-full-pool-20260926.json"
inputs = {
    "baseline-config.json": source / "config.json",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "measurements.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
}
names = set(baseline["foundation"] + baseline["mordred"])
names.add(baseline["weights"]["classifier"])
for fold in baseline["folds"]:
    names.add(fold["oof"])
    paths = fold["endpoints"]["CYP3A4_is_TDI"]
    names.update(paths[key] for key in ("model", "columns", "imputer", "pca"))
    for parent in fold["auxiliary"]:
        names.add(parent["summary"])
        names.update(parent["predictions"])
inputs.update({f"old-{name}": source / "inputs" / name for name in sorted(names)})
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
    "Extend the exact submitted A4 historical TabPFN stack to all UMAP held-fold molecules",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
completed = subprocess.run(
    command + ["--", "python", "scripts/classification_a4_historical_full_pool.py"],
    cwd=root,
    check=False,
)
if completed.returncode:
    raise SystemExit(f"Historical A4 full-pool run failed with exit code {completed.returncode}")
