"""Record the intrinsic CYP continuous architecture screen."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
CONFIG = ROOT / "configs/external-training/task-focus-cyp-architecture-v1.json"
command = [
    "pixi", "run", "cyp", "run", "--track", "shared",
    "--purpose", "Intrinsic architecture screen for the fixed 14-task CYP continuous model",
    "--config", str(CONFIG),
    "--input", f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
    "--input", f"observations.parquet={CORPUS / 'observations-endpoint-only.parquet'}",
    "--input", f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
    "--", "python", "scripts/tune_external_task_architecture.py",
]
subprocess.run(command, cwd=ROOT, check=True)
