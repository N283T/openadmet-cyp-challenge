"""Launch BioTransformer target-association and importance diagnostics."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
classification = root / "artifacts/runs/classification"
inputs = [
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    classification
    / "20260915T065421221917Z-7d41918f/outputs/product-activity-features.parquet",
    classification
    / "20260915T070040196804Z-517e5bfc/outputs/product-embeddings.parquet",
    classification
    / "20260915T070040196804Z-517e5bfc/outputs/product-embedding-scalars.parquet",
]
for path in inputs:
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
    "Diagnose TDI target associations and held-out importance of metabolite features",
    "--config",
    str(root / "configs/classification-round/biotransformer-signal-diagnostics-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command
    + ["--", "python", "scripts/classification_biotransformer_signal_diagnostics.py"],
    cwd=root,
    check=True,
)
