"""Launch the matched regression dual-split feature audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit feature-candidate ordering across parent-retained series and scaffold splits",
    "--config",
    str(root / "configs/regression/dual-split-feature-audit-v1.json"),
    "--input",
    f"primary-metrics.csv={runs / '20260918T012518053809Z-5c46e72d/outputs/metrics.csv'}",
    "--input",
    f"secondary-metrics.csv={runs / '20260918T021832961989Z-c89b25fa/outputs/metrics.csv'}",
    "--",
    "python",
    "scripts/regression_dual_split_feature_audit.py",
]
subprocess.run(command, cwd=root, check=True)
