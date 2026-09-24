"""Record D6 call-count audit for control and conservative blend."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
config = root / "configs/analysis/d6-call-count-audit-v1.json"
inputs = {
    "blend-oof.parquet": runs
    / "classification/20260924T043241076761Z-7352d993/outputs/all-pool-oof.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "blind-scores.parquet": runs
    / "classification/20260923T173655387004Z-97d71dce/outputs/decision-test.parquet",
    "base-submission.parquet": runs
    / "classification/20260915T125939829893Z-1d921ffc/outputs/predictions.parquet",
    "aliases.parquet": runs / "shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
}
declared = json.loads(config.read_text())
expected = {
    declared["oof"],
    declared["assignments"],
    declared["blind_scores"],
    declared["base_submission"],
    declared["aliases"],
}
if set(inputs) != expected:
    raise ValueError("Config and input aliases differ")
for path in [config, *inputs.values()]:
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
    "Audit D6 model and positive-call count jointly across complete-pool OOF and blind scores",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_call_count_audit.py"],
    cwd=root,
    check=True,
)
