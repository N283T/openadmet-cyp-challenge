"""Extend the frozen prior global pool with completed 2D6-focused experiments."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs/regression"
parser = argparse.ArgumentParser()
parser.add_argument("--tag", required=True)
args = parser.parse_args()
assert re.fullmatch(r"[a-z0-9_]+", args.tag)
prior = b / "20260912T094021204217Z-84c28266"
config = json.loads((prior / "config.json").read_text())
sources = config["sources"]
new = []
exclude_runs = {
    "20260912T101320386840Z-dcf09b87": "Unconverged alpha solver; superseded by320cdd53"
}
allowed = {
    "d6_kernel.py",
    "d6_interval_kernel.py",
    "d6_efficacy_readout.py",
    "d6_interval_width.py",
    "d6_specialist_fit.py",
    "d6_specialist_seed_summary.py",
    "d6_joint_external.py",
    "d6_joint_external_fit.py",
}
for run in sorted(b.iterdir()):
    if run.name < "20260912T100000" or run.name in exclude_runs:
        continue
    if not (run / "result.json").exists():
        continue
    if json.loads((run / "result.json").read_text()).get("status") != "completed":
        continue
    meta = json.loads((run / "manifest.json").read_text())
    script = next((Path(s).name for s in meta.get("command", []) if s.endswith(".py")), "")
    cfg = json.loads((run / "config.json").read_text()) if (run / "config.json").exists() else {}
    eligible = script in allowed or (
        script == "focus_single_readout.py"
        and str(cfg.get("method", "")).startswith("d6_specialist")
    )
    if not eligible:
        continue
    for p in sorted((run / "outputs").glob("*oof*.parquet")):
        if p.name.startswith("canonical"):
            continue
        new.append(dict(path=str(p), run=run.name, script=script, purpose=meta["purpose"]))
assert len(new) >= 10, len(new)
config["sources"] = sources + new
config["explicit_excluded_runs"] = exclude_runs
config["policy"] += (
    "; frozen prior pool plus completed 2D6-focused round; exclude unconverged solver"
)
p = root / f"configs/global-caruana/{args.tag}.json"
with p.open("x") as f:
    json.dump(config, f, indent=2)
scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"
reference = b / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet"


def launch(script, inputs, purpose):
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        purpose,
        "--config",
        str(p),
    ]
    for f in inputs:
        cmd += ["--input", str(f)]
    result = subprocess.run(
        cmd + ["--", "python", script], cwd=root, text=True, capture_output=True
    )
    print(result.stdout, result.stderr, flush=True)
    result.check_returncode()
    return re.search(r"run_id=(\S+)", result.stdout)[1]


run = launch(
    "scripts/global_caruana.py",
    [reference, scorer] + [Path(s["path"]) for s in config["sources"]],
    "Expand frozen global candidate pool with audited 2D6 focused experiments",
)
d = b / run / "outputs"
audit = launch(
    "scripts/audit_global_caruana.py",
    [d / "aligned-pool.parquet", d / "global-oof.parquet", d / "weights.json", scorer],
    "Independently reconstruct expanded 2D6 Caruana weights and scores",
)
with (root / f"artifacts/launches/{args.tag}_manifest.json").open("x") as f:
    json.dump(dict(selection_run=run, audit_run=audit, new_sources=new), f, indent=2)
