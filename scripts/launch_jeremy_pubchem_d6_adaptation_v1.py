"""Launch Jeremy-initialized PubChem-only CYP2D6 domain adaptation."""

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
CONNECTIVITY = (
    ROOT
    / "artifacts/runs/shared/20260911T052745709989Z-8792f422/outputs"
    / "challenge-connectivity.parquet"
)
OWNED_RUN = ROOT / "artifacts/runs/regression/20260916T143736207347Z-cc6eb40b/outputs"
JEREMY_RUN = ROOT / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
JEREMY_CHECKPOINT = Path(
    "/home/nagaet/ghq/github.com/jeremycheminf/openadmet_scripts/"
    "CYP_Challenge/checkpoints/chemprop_medium.pt"
)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test Jeremy-initialized PubChem AID 1851-only domain adaptation for CYP2D6",
    "--config",
    str(ROOT / "configs/external-training/jeremy-pubchem-d6-adaptation-v1.json"),
    "--input",
    f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
    "--input",
    f"observations-compact-all.parquet={CORPUS / 'observations-compact-all.parquet'}",
    "--input",
    f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
    "--input",
    f"challenge-connectivity.parquet={CONNECTIVITY}",
    "--input",
    f"measurements.parquet={CHALLENGE / 'measurements.parquet'}",
    "--input",
    f"parent-assignments.parquet={ASSIGNMENTS}",
    "--input",
    "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
    "--input",
    f"jeremy-checkpoint.pt={JEREMY_CHECKPOINT}",
    "--input",
    "owned-oof.parquet="
    + str(OWNED_RUN / "trial-03-joint_all_no_binary/external-oof.parquet"),
    "--input",
    f"jeremy-oof.parquet={JEREMY_RUN / 'external-oof.parquet'}",
    "--",
    "python",
    "scripts/train_embedding_distillation.py",
]
subprocess.run(command, cwd=ROOT, check=True)
