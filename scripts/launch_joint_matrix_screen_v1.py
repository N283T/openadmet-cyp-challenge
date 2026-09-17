"""Launch the project-owned joint sparse target-matrix screen."""

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
CHALLENGE_CONNECTIVITY = (
    ROOT
    / "artifacts/runs/shared/20260911T052745709989Z-8792f422/outputs"
    / "challenge-connectivity.parquet"
)
CONFIG = ROOT / "configs/external-training/joint-matrix-screen-v1.json"

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Train project-owned joint sparse-matrix encoders and audit target-specific "
    "raw LightGBM transfer",
    "--config",
    str(CONFIG),
    "--input",
    f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
    "--input",
    f"observations-assay.parquet={CORPUS / 'observations-assay.parquet'}",
    "--input",
    f"observations-compact-all.parquet={CORPUS / 'observations-compact-all.parquet'}",
    "--input",
    f"observations-compact-cyp.parquet={CORPUS / 'observations-compact-cyp.parquet'}",
    "--input",
    f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
    "--input",
    f"challenge-connectivity.parquet={CHALLENGE_CONNECTIVITY}",
    "--input",
    f"measurements.parquet={CHALLENGE / 'measurements.parquet'}",
    "--input",
    f"parent-assignments.parquet={ASSIGNMENTS}",
    "--input",
    "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
    "--input",
    "previous-owned.parquet="
    + str(
        ROOT
        / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
        / "external-oof.parquet"
    ),
    "--input",
    "jeremy-oof.parquet="
    + str(
        ROOT
        / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
        / "external-oof.parquet"
    ),
    "--",
    "python",
    "scripts/train_external_joint_matrix.py",
]
subprocess.run(command, cwd=ROOT, check=True)
