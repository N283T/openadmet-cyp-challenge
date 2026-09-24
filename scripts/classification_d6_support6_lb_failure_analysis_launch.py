"""Launch the recorded D6 support6 leaderboard-failure diagnosis."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
support_run = runs / "20260918T085333214169Z-892612f4"
consensus_run = runs / "20260918T085712373599Z-8829d528"
pre_lb_run = root / "artifacts/api/classification/20260918T094508388319Z-af6aa011"
post_lb_run = root / "artifacts/api/classification/20260918T101216302522Z-4a8b944c"


def downloaded(run: Path) -> Path:
    path = Path(json.loads((run / "response.json").read_text())["result"])
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


inputs = {
    "consensus-test.parquet": consensus_run / "outputs/test.parquet",
    "consensus-oof.parquet": consensus_run / "outputs/oof.parquet",
    "member-test.parquet": consensus_run / "inputs/split-bag-test.parquet",
    "member-oof.parquet": consensus_run / "inputs/split-bag-oof.parquet",
    "oof-support.parquet": consensus_run / "inputs/random-oof-support.parquet",
    "test-support.parquet": support_run / "outputs/test-support.parquet",
    "compounds.parquet": support_run / "inputs/compounds.parquet",
    "measurements.parquet": support_run / "inputs/measurements.parquet",
    "test-mapping.parquet": support_run / "inputs/test-mapping.parquet",
    "pre-d6-leaderboard.csv": downloaded(pre_lb_run),
    "post-d6-leaderboard.csv": downloaded(post_lb_run),
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
    "Diagnose the failed support6 CYP2D6 leaderboard probe at swapped-compound grain",
    "--config",
    str(root / "configs/analysis/d6-support6-lb-failure-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command
    + ["--", "python", "scripts/classification_d6_support6_lb_failure_analysis.py"],
    cwd=root,
    check=True,
)
