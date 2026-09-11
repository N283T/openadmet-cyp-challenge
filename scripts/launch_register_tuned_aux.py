"""Record validation of adopted auxiliary assets and publish their explicit registry."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--revision", default="")
args = parser.parse_args()
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = {
    "selection_run": "20260911T120540987924Z-3f5823d1",
    "split_run": "20260911T100134004950Z-c7a3f2bf",
    "compounds": add(
        root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
    ),
    "assignments": add(
        root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet"
    ),
    "models": [],
    "selected": add(
        root / "artifacts/runs/shared/20260911T120540987924Z-3f5823d1/outputs/selected.json"
    ),
}
for log in ["tuned-aux-20260911.log", "tuned-aux-remaining-20260911.log"]:
    for run in re.findall(r"run_id=(\S+)", (root / "artifacts/launches" / log).read_text()):
        out = root / "artifacts/runs/shared" / run / "outputs"
        cfg["models"].append(
            {
                "run": run,
                "summary": add(out / "summary.json"),
                "membership": add(out / "membership.parquet"),
                "predictions": [add(p) for p in sorted(out.glob("predictions-*.parquet"))],
            }
        )
p = root / f"configs/parent-cv/adopt-tuned-aux{args.revision}.json"
with p.open("x") as stream:
    json.dump(cfg, stream, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Validate and register all ten adopted Rustuna auxiliary models",
    "--config",
    str(p),
]
for path in paths:
    cmd += ["--input", str(path)]
result = subprocess.run(
    cmd + ["--", "python", "scripts/register_tuned_aux.py"],
    check=False,
    cwd=root,
    capture_output=True,
    text=True,
)
run = re.findall(r"run_id=(\S+)", result.stdout)[-1]
if result.returncode:
    raise RuntimeError(f"Registration failed; inspect artifacts/runs/shared/{run}/console.log")
registry = root / "artifacts/runs/shared" / run / "outputs/registry.json"
with (root / "configs/regression/current-log2fc.json").open("x") as stream:
    stream.write(registry.read_text())
print("run_id=" + run)
