"""Record the Jeremy-medium overlap and representation diagnostic."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
JEREMY = ROOT / "artifacts/runs/shared/20260912T055819206227Z-c060d197/outputs"
OBJECTIVE = ROOT / "artifacts/runs/regression/20260916T120842454770Z-b5a003ef/outputs"
SOURCE = ROOT / "artifacts/runs/regression/20260916T123149468017Z-68f542a6/outputs"
FULL_OOF = ROOT / "artifacts/runs/regression/20260916T125409299702Z-7e362fe6/outputs"
PREVIOUS = ROOT / "artifacts/runs/shared/20260916T085412850783Z-b1f5762b/outputs"
CONFIG = ROOT / "configs/owned-encoder-v2/jeremy_medium_gap_diagnostic_v1.json"

inputs = {
    "parent-assignments.parquet": ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet",
    "challenge-overlap-observations.parquet": CORPUS / "challenge-overlap-observations.parquet",
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
    "jeremy-oof.parquet": ROOT
    / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs/external-oof.parquet",
    "owned-1a2.parquet": FULL_OOF / "source_sqrt_no_tox21/external-oof.parquet",
    "owned-2c9.parquet": FULL_OOF / "source_sqrt_no_polaris/external-oof.parquet",
    "owned-2d6.parquet": FULL_OOF / "objective_sqrt_no_classification/external-oof.parquet",
    "owned-3a4.parquet": FULL_OOF / "objective_proportional_mixed/external-oof.parquet",
    "jeremy-embeddings.parquet": JEREMY / "external-embeddings.parquet",
    "owned-sqrt-embeddings.parquet": OBJECTIVE / "trial-02-sqrt_mixed/external-embeddings.parquet",
    "owned-previous-embeddings.parquet": PREVIOUS / "external-embeddings.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Diagnose Jeremy-medium representation gap and possible public-source overlap",
    "--config",
    str(CONFIG),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/analyze_jeremy_medium_gap.py"])
subprocess.run(command, cwd=ROOT, check=True)
