"""Declare canonical lineage, raw source and current official leaderboard evidence."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = []
cfg = {}
for key, name in [
    ("aliases", "aliases"),
    ("standard", "standardization"),
    ("labels", "measurements"),
    ("conflicts", "measurement_conflicts"),
    ("compounds", "compounds"),
]:
    p = base / (name + ".parquet")
    cfg[key] = f"{len(inputs):03d}-{p.name}"
    inputs.append(p)
cfg["raw"] = []
for p in sorted((root / "data/raw").glob("*TRAIN*.csv")):
    if "single-concentration" in p.name:
        continue
    cfg["raw"].append(f"{len(inputs):03d}-{p.name}")
    inputs.append(p)
cfg["leaderboards"] = {}
for ep, run in [
    ("MA", "20260912T125633423600Z-738ff4c9"),
    ("CYP2D6", "20260912T125633429611Z-e4ea4a7a"),
    ("CYP3A4", "20260912T125633453349Z-64e5e66c"),
]:
    p = next((root / "artifacts/api/classification" / run / "downloads").glob("*/*.csv"))
    cfg["leaderboards"][ep] = f"{len(inputs):03d}-{p.name}"
    inputs.append(p)
for p in [
    root / "artifacts/metric-reference/config.py",
    root / "artifacts/metric-reference/evaluate_predictions.py",
    root / "artifacts/api/20260912T113927495541Z-f38c3640/leaderboards.py",
    root / "references/official/cyp/challenge/20260907T135802Z/snapshot.md",
    root / "docs/submission-baseline.md",
]:
    inputs.append(p)
path = root / "configs/classification/raw-provenance.json"
path.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Reconcile supplied raw TDI booleans to canonical labels and verify "
    "current official classification baseline numbers",
    "--config",
    str(path),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_label_provenance.py"], check=True, cwd=root
)
