"""Declare paired fold-zero predictions and auxiliary lineage."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--readout-run", required=True)
parser.add_argument("--search-run", required=True)
args = parser.parse_args()
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = {
    "original": add(
        root
        / "artifacts/runs/regression/20260911T105206572152Z-e998787e/outputs/parent-oof.parquet"
    ),
    "tuned": add(root / f"artifacts/runs/regression/{args.readout_run}/outputs/parent-oof.parquet"),
    "scorer": add(root / "artifacts/metric-reference/custom_scoring_functions.py"),
    "targets": add(
        root / "artifacts/runs/shared/20260910T100337258849Z-927eac1c/outputs/targets.parquet"
    ),
    "assignments": add(
        root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet"
    ),
    "aux": {},
    "database": {"host": str(root / ".runtime/socket"), "port": 5434, "dbname": "cyp_challenge"},
}
search = root / "artifacts/runs/shared" / args.search_run / "outputs"
winner = json.loads((search / "selected.json").read_text())
cfg["trials"] = add(search / "trials.json")
cfg["selected_trial"] = winner["trial"]
cfg["histories"] = {
    "original": add(search / "trial-000/history.jsonl"),
    "selected": add(search / f"trial-{winner['trial']:03d}/history.jsonl"),
}
for method, log in [
    ("original", "parent-cv-aux-20260911.log"),
    ("tuned", "tuned-aux-20260911.log"),
]:
    cfg["aux"][method] = []
    for run in re.findall(r"run_id=(\S+)", (root / "artifacts/launches" / log).read_text()):
        out = root / "artifacts/runs/shared" / run / "outputs"
        meta = json.loads((out / "summary.json").read_text())
        if 0 in meta["config"]["excluded_folds"]:
            cfg["aux"][method].append([add(p) for p in sorted(out.glob("predictions-*.parquet"))])
    assert len(cfg["aux"][method]) == 4
p = root / "configs/parent-cv/tuned-audit.json"
with p.open("x") as stream:
    json.dump(cfg, stream, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare Rustuna-selected auxiliary versus original on untouched outer fold zero",
    "--config",
    str(p),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(cmd + ["--", "python", "scripts/audit_tuned_log2fc.py"], check=True, cwd=root)
