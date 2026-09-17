"""Audit first-wave v2 corpus variants against v1 owned and Jeremy encoders."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/owned-encoder-v2/data-screen-audit.json"
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
RUNS = {
    "v2_assay": "20260916T082842598129Z-546bf741",
    "v2_compact_all": "20260916T083209317276Z-d8372328",
    "v2_compact_cyp": "20260916T083542263535Z-681f18b0",
    "v2_endpoint_only": "20260916T083915098714Z-00da8699",
    "v1_owned_w600": "20260916T075207221915Z-2bb17166",
    "jeremy_raw600_lgbm": "20260912T055849013223Z-f16a5dfd",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Matched all-target raw LightGBM audit of v2 corpus variants",
    "--config",
    str(CONFIG),
    "--input",
    f"custom_scoring_functions.py={SCORER}",
]
for name, run in RUNS.items():
    path = ROOT / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet"
    command.extend(["--input", f"{name}.parquet={path}"])
command.extend(["--", "python", "scripts/audit_atomic_readout_grid.py"])
subprocess.run(command, cwd=ROOT, check=True)
