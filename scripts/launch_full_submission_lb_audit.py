"""Record official leaderboard comparison for the submitted full-data ensemble."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg = dict(submission_id="20260912T074802231782Z-eaafa860", snapshots=[])
inputs = []
for ep, run, oof in [
    ("MA", "20260912T081048806309Z-6fe81b63", 0.5964141092282933),
    ("CYP1A2", "20260912T081116173280Z-1750752a", 0.683335),
    ("CYP2C9", "20260912T081116169013Z-c3fe0073", 0.463781),
    ("CYP2D6", "20260912T081116179549Z-5b722cbb", 0.834741),
    ("CYP3A4", "20260912T081116193953Z-d9b548b8", 0.403800),
]:
    paths = list((root / f"artifacts/api/regression/{run}/downloads").rglob("*.csv"))
    assert len(paths) == 1
    cfg["snapshots"].append(
        dict(endpoint=ep, snapshot=run, oof=oof, file=f"{len(inputs):03d}-{paths[0].name}")
    )
    inputs.append(paths[0])
config = root / "configs/submission-full/leaderboard-first-result.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Match submitted timestamp and compare full-data ensemble LB with official TabICL",
    "--config",
    str(config),
]
for path in inputs:
    cmd += ["--input", str(path)]
subprocess.run(cmd + ["--", "python", "scripts/audit_full_submission_lb.py"], cwd=root, check=True)
