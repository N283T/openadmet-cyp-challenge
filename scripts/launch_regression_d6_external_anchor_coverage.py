"""Record external CYP2D6 measurement coverage near challenge hit anchors."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/shared"
inputs = {
    "parent-assignments.parquet": (
        root / "artifacts/runs/regression"
        / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "challenge-compounds.parquet": (
        runs / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
    ),
    "external-compounds.parquet": (
        runs / "20260916T082302554619Z-656c63e6/outputs/external-compounds.parquet"
    ),
    "external-observations.parquet": (
        runs / "20260916T082302554619Z-656c63e6/outputs/observations-assay.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Audit independent external D6 assay coverage near other-CYP hit anchors",
    "--config", str(root / "configs/regression/d6-external-anchor-coverage-v3.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_external_anchor_coverage.py"],
    cwd=root,
    check=True,
)
