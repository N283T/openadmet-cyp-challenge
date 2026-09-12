"""Declare completed OOF artifacts before global ensemble selection."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--tag", default="all_saved_v1")
args = parser.parse_args()
assert re.fullmatch(r"[a-z0-9_]+", args.tag)
base = root / "artifacts/runs/regression"
allowed_audits = {
    "audit_readouts.py",
    "audit_tabicl.py",
    "audit_tabicl_parent.py",
    "audit_parent_cv.py",
    "audit_foundation_replacement.py",
}
excluded = (
    "audit_",
    "calibrat",
    "compare_",
    "fixed_",
    "save_",
    "summarize_",
    "diagnosis",
    "distribution",
    "context_",
    "infer_blend",
    "global_caruana",
)
files = []
for run in sorted(base.iterdir()):
    result = run / "result.json"
    manifest = run / "manifest.json"
    if not result.exists() or not manifest.exists():
        continue
    if json.loads(result.read_text()).get("status") != "completed":
        continue
    meta = json.loads(manifest.read_text())
    script = next((Path(s).name for s in meta.get("command", []) if s.endswith(".py")), "")
    if not script or (script not in allowed_audits and any(s in script for s in excluded)):
        continue
    for p in sorted((run / "outputs").glob("*oof*.parquet")):
        if any(s in p.name for s in ["progress", "aligned", "calibrat"]):
            continue
        files.append(dict(path=str(p), run=run.name, script=script, purpose=meta["purpose"]))
reference = base / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet"
scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"
config = dict(
    steps=50,
    seed=42,
    sources=files,
    policy=(
        "Completed compatible OOF artifacts, including prior ensembles; no score prefilter; "
        "exact endpoint keys and labels; duplicate vectors collapsed"
    ),
    crossfit_interpretation=(
        "Weight cross-fit only; base models are not retrained in nested meta folds"
    ),
)
p = root / f"configs/global-caruana/{args.tag}.json"
with p.open("x") as f:
    json.dump(config, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Global Caruana selection across all compatible saved model-family OOF predictions",
    "--config",
    str(p),
]
for f in [reference, scorer] + [Path(s["path"]) for s in files]:
    cmd += ["--input", str(f)]
print(f"Declared {len(files)} candidate artifacts", flush=True)
r = subprocess.run(
    cmd + ["--", "python", "scripts/global_caruana.py"], cwd=root, text=True, capture_output=True
)
print(r.stdout, r.stderr, flush=True)
r.check_returncode()
run_id = re.search(r"run_id=(\S+)", r.stdout)[1]
paths = [
    base / run_id / "outputs" / name
    for name in ["aligned-pool.parquet", "global-oof.parquet", "weights.json"]
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Independently reconstruct global Caruana weights and official OOF metrics",
    "--config",
    str(p),
]
for f in paths + [scorer]:
    cmd += ["--input", str(f)]
r = subprocess.run(
    cmd + ["--", "python", "scripts/audit_global_caruana.py"],
    cwd=root,
    text=True,
    capture_output=True,
)
print(r.stdout, r.stderr, flush=True)
r.check_returncode()
