"""Record the three-submission fold-excluded primary-screen proxy audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/analysis/d6-threeway-primary-proxy-v1.json"
inputs = {
    "threeway-test.parquet": root
    / "artifacts/runs/classification/20260923T190811462293Z-4f4a92f2/outputs/threeway-test.parquet",
    "threeway-oof.parquet": root
    / "artifacts/runs/classification/20260923T190811462293Z-4f4a92f2/outputs/threeway-oof.parquet",
    "all-nested-log2fc.parquet": root
    / "artifacts/runs/shared/20260912T082554377951Z-d272774c/outputs/all-nested-log2fc.parquet",
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
    "Test independent fold-excluded primary-screen proxy on September D6 swaps",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_threeway_primary_proxy.py"],
    cwd=root,
    check=True,
)
