"""Launch the recorded Morgan attribution summary."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = root / "artifacts/runs/classification/20260913T124832303603Z-67007a7a/outputs"
paths = [source / "bit-attributions.csv", source / "ablation-metrics.csv"]
assert all(path.is_file() for path in paths)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Extract fold-stable cross-model Morgan-bit explanations from attribution outputs",
    "--config",
    str(root / "configs/classification-round/fp-attribution-summary-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_fp_attribution_summary.py"],
    cwd=root,
    check=True,
)
