"""Record complete-pool OOF for the frozen MACE and SMARTCyp AP atoms."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
classification = root / "artifacts/runs/classification"
config = root / "configs/classification/d6-ap-mechanism-atoms-full-pool-v1.json"
inputs = {
    "compounds.parquet": shared / "compounds.parquet",
    "measurements.parquet": shared / "measurements.parquet",
    "test-mapping.parquet": shared / "test_mapping.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "mace-off-embeddings.parquet": classification
    / "20260914T050016310677Z-e80256e2/outputs/mace-off-embeddings.parquet",
    "smartcyp-features.parquet": classification
    / "20260914T044537363003Z-5e7f9e02/outputs/smartcyp-features.parquet",
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
    "Reproduce frozen MACE and SMARTCyp AP atoms on complete held-fold pools",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_embedding_readout.py"],
    cwd=root,
    check=True,
)
