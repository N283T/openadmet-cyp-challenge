"""Record a UNIQUE-style D6 applicability-domain audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "parent-assignments.parquet": (
        runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "scaffold-assignments.parquet": (
        runs / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet"
    ),
    "measurements.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
    "compounds.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
    ),
    "test-mapping.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"
    ),
    "medium.parquet": (
        runs / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet"
    ),
    "parent-oof.parquet": (
        runs / "regression/20260926T013148480112Z-72cf4682/outputs/aligned-oof.parquet"
    ),
    "scaffold-oof.parquet": (
        runs / "regression/20260918T021832961989Z-c89b25fa/outputs/jeremy_anchor_all-oof.parquet"
    ),
    "parent-members.parquet": (
        root / "artifacts/tracked/d6-depth-atomic-20260925/ensemble-oof.parquet"
    ),
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Measure D6 applicability-domain indicators against fold-held error and blind chemistry",
    "--config",
    str(root / "configs/regression/d6-applicability-domain-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_applicability_domain.py"],
    cwd=root,
    check=True,
)
