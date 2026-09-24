"""Launch a recorded classification candidate preflight from explicit inputs."""

import argparse
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--config", type=Path, required=True)
parser.add_argument("--decision-test", type=Path, required=True)
parser.add_argument("--decision-oof", type=Path, required=True)
parser.add_argument("--member-test", type=Path, required=True)
parser.add_argument("--member-oof", type=Path, required=True)
parser.add_argument("--compounds", type=Path, required=True)
parser.add_argument("--measurements", type=Path, required=True)
parser.add_argument(
    "--purpose",
    default="Run fail-closed classification candidate preflight",
)
args = parser.parse_args()

config = args.config.resolve()
inputs = {
    "decision-test.parquet": args.decision_test.resolve(),
    "decision-oof.parquet": args.decision_oof.resolve(),
    "member-test.parquet": args.member_test.resolve(),
    "member-oof.parquet": args.member_oof.resolve(),
    "compounds.parquet": args.compounds.resolve(),
    "measurements.parquet": args.measurements.resolve(),
}
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
    args.purpose,
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_candidate_preflight.py"],
    cwd=root,
    check=True,
)
