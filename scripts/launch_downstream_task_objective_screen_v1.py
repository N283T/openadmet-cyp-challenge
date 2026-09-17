"""Record the downstream-guided external task-objective screen."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
CONFIG = ROOT / "configs/external-training/downstream-task-objective-screen-v1.json"
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Downstream-guided target-specific screen of external task objectives",
    "--config", str(CONFIG),
    "--input", f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
    "--input", f"observations.parquet={CORPUS / 'observations-endpoint-only.parquet'}",
    "--input", f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
    "--input", f"measurements.parquet={CHALLENGE / 'measurements.parquet'}",
    "--input", f"parent-assignments.parquet={ASSIGNMENTS}",
    "--input", "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
    "--input", "previous-owned.parquet="
    + str(
        ROOT
        / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
        / "external-oof.parquet"
    ),
    "--input", "jeremy-oof.parquet="
    + str(
        ROOT
        / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
        / "external-oof.parquet"
    ),
    "--", "python", "scripts/screen_downstream_task_objectives.py",
]
subprocess.run(command, cwd=ROOT, check=True)
