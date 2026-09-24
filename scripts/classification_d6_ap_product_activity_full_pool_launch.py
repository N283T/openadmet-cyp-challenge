"""Record complete-pool OOF for the frozen product activity AP atom."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
config = root / "configs/classification/d6-ap-product-activity-full-pool-v1.json"
inputs = {
    "compounds.parquet": shared / "compounds.parquet",
    "measurements.parquet": shared / "measurements.parquet",
    "test-mapping.parquet": shared / "test_mapping.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "product-activity-features.parquet": (
        root
        / "artifacts/runs/classification/20260915T065421221917Z-7d41918f/outputs/"
        "product-activity-features.parquet"
    ),
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
    "Reproduce frozen product activity AP atom on complete held-fold pools",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_biotransformer_activity_readout.py"],
    cwd=root,
    check=True,
)
