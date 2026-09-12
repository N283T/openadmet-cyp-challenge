"""Declare exact source tables for follow-up diagnostics, with no submission."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    p = root / path
    assert p.is_file()
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


run = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
base = json.loads((run / "config.json").read_text())
cfg = dict(
    oof=add(run / "outputs/ensemble-oof.parquet"),
    old=add("artifacts/runs/shared/20260910T123856868493Z-243a1e23/outputs/oof.parquet"),
    oracle=add(
        "artifacts/runs/classification/20260912T125855315676Z-6559a18f/outputs/oracle_thresholds.csv"
    ),
    cuts=add(run / "outputs/thresholds.csv"),
    ensemble_config=add(run / "config.json"),
    sources=[
        dict(prefix=s["prefix"], test=add(run / "inputs" / s["test"])) for s in base["sources"]
    ],
    selected={"CYP2D6_is_TDI": "all_plain_pfn_graph_mean", "CYP3A4_is_TDI": "lgbm_balanced_mean"},
)
config = root / "configs/classification-round/diagnostic-followup-v1.json"
assert not config.exists()
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Paired OOF uncertainty and fold threshold deployment diagnostic without refitting",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/tdi_diagnostic_followup.py"], cwd=root, check=True)
