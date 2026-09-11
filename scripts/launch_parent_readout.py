"""Launch complete-feature parent-retention comparison after nested auxiliary fits."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
parser = argparse.ArgumentParser()
parser.add_argument("--fold", type=int, action="append")
parser.add_argument(
    "--aux-log", help="Explicit historical launch log instead of the current registry"
)
parser.add_argument("--aux-registry", default="configs/regression/current-log2fc.json")
parser.add_argument("--revision", default="")
args = parser.parse_args()
parents = []
if args.aux_log:
    runs = re.findall(r"run_id=(\S+)", (root / args.aux_log).read_text())
else:
    runs = [m["run"] for m in json.loads((root / args.aux_registry).read_text())["models"]]
for run in runs:
    out = shared / run / "outputs"
    if not (out / "summary.json").exists():
        continue
    meta = json.loads((out / "summary.json").read_text())
    parents.append((out, meta))
for fold in args.fold if args.fold is not None else range(5):
    assert sum(fold in meta["config"]["excluded_folds"] for _, meta in parents) == 4
    paths = []

    def add(p, paths=paths):
        paths.append(p)
        return f"{len(paths) - 1:03d}-{p.name}"

    add(
        root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet"
    )
    add(shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet")
    for run, pattern in [
        ("20260910T012630927160Z-e9576045", "chemeleon-*.parquet"),
        ("20260909T132918829974Z-bc81af58", "mordred-*.parquet"),
    ]:
        for p in sorted((shared / run / "outputs").glob(pattern)):
            add(p)
    cfg = dict(
        fold=fold,
        split_run="20260911T104703655337Z-a3923255",
        auxiliary=[],
        model="LightGBM L1 300",
        seed=42,
        features="CheMeleon PCA256 + nested CheMeleon log2fc4 + Mordred; descriptor ablation",
    )
    for out, meta in parents:
        if fold in meta["config"]["excluded_folds"]:
            cfg["auxiliary"].append(
                dict(
                    summary=add(out / "summary.json"),
                    membership=add(out / "membership.parquet"),
                    predictions=[add(p) for p in sorted(out.glob("predictions-*.parquet"))],
                )
            )
    path = root / f"configs/parent-cv/readout-{fold}{args.revision}.json"
    with path.open("x") as f:
        json.dump(cfg, f, indent=2)
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Paired parent-retained/excluded full-feature LightGBM comparison fold {fold}",
        "--config",
        str(path),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/parent_readout.py"], check=True, cwd=root)
