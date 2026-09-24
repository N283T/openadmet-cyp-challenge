"""Launch the random-OOF to blind-test-support post-stratification audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "compounds.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "tabpfn-random-oof.parquet": runs
    / "classification/20260918T084809611158Z-8b586046/outputs/oof.parquet",
    "tree-random-oof.parquet": runs
    / "classification/20260918T085057522996Z-3f5b3fc3/outputs/oof.parquet",
    "deployed-oof.parquet": runs
    / "classification/20260918T063800192024Z-c1058115/outputs/oof.parquet",
}
for path in inputs.values():
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Post-stratify random-fold CYP2D6 models to the blind-test support distribution",
    "--config",
    str(root / "configs/classification/d6-random-support-reweight-audit-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_random_support_reweight_audit.py"],
    cwd=root,
    check=True,
)
