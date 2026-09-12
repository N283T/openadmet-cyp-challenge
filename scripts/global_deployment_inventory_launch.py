"""Record metadata-only deployment inventory; never fit or predict."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/regression"
global_run = base / "20260912T094021204217Z-84c28266"
catalog = json.loads((global_run / "outputs/catalog.json").read_text())
weights = json.loads((global_run / "outputs/weights.json").read_text())
selected = {
    k
    for w in weights
    if w["method"] == "caruana_full_fit" and "2D6" not in w["endpoint"]
    for k in w["weights"]
}
runs = {r["source"]["run"] for r in catalog if r["column"] in selected}
# Resolve nested fixed means and the earlier audit-only LightGBM baseline chain.
pending = list(runs)
while pending:
    name = pending.pop()
    cfg = json.loads((base / name / "config.json").read_text())
    children = [
        c["run_id"] for c in cfg.get("components", []) if isinstance(c, dict) and "run_id" in c
    ]
    if name == "20260911T111952347762Z-e9f6c457":
        children += ["20260911T110118749415Z-5167f1a0"]
    if name == "20260911T110118749415Z-5167f1a0":
        meta = json.loads((base / name / "manifest.json").read_text())
        for i in meta["inputs"]:
            p = Path(i["original"])
            if p.name.endswith("oof.parquet") and p.parent.name == "outputs":
                children.append(p.parent.parent.name)
    for child in children:
        if child not in runs:
            runs.add(child)
            pending.append(child)
inputs = [global_run / "outputs/catalog.json", global_run / "outputs/weights.json"]
records = []
for name in sorted(runs):
    files = {}
    for kind in ["manifest", "config", "result"]:
        p = base / name / (kind + ".json")
        files[kind] = f"{len(inputs):03d}-{p.name}"
        inputs.append(p)
    records.append(dict(run=name, files=files))
config = root / "configs/scp-calibration/global-deployment-inventory.json"
config.write_text(
    json.dumps(
        dict(
            global_run=global_run.name,
            records=records,
            policy="Read-only metadata inventory; no fit or inference; "
            "global126 full-OOF weight support for1A2/2C9/3A4",
        ),
        indent=2,
    )
)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Inventory selected global126 model and nested aggregate deployment support; "
    "no training or inference",
    "--config",
    str(config),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/global_deployment_inventory.py"], check=True, cwd=root
)
