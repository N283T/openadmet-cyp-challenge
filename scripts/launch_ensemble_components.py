"""Launch matched three-component OOF experiments with explicit snapshots."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
registry = json.loads((root / "configs/regression/current-parent-readout.json").read_text())
for item in registry["models"]:
    parent = root / "artifacts/runs/regression" / item["run"]
    cfg = json.loads((parent / "config.json").read_text())
    paths = []

    def add(p, paths=paths):
        paths.append(p)
        return f"{len(paths) - 1:03d}-{p.name}"

    def remap(obj, parent=parent):
        if isinstance(obj, dict):
            return {k: remap(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [remap(v) for v in obj]
        if isinstance(obj, str) and (parent / "inputs" / obj).is_file():
            return add(parent / "inputs" / obj)
        return obj

    cfg = remap(cfg)
    for p in sorted((parent / "inputs").iterdir()):
        if p not in paths:
            add(p)
    add(root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet")
    cfg["weights"] = add(Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt")
    cfg["components"] = [
        "Morgan2048 radius2 + Mordred -> LGBM L1 300",
        "frozen CheMeleon2048 -> TabPFN2.6",
        "nested pred_log2fc4 -> TabPFN2.6",
    ]
    p = root / f"configs/parent-cv/ensemble-components-v2-{item['fold']}.json"
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
        f"Fit three separate ensemble components on matched fold {item['fold']}",
        "--config",
        str(p),
    ]
    for f in paths:
        cmd += ["--input", str(f)]
    subprocess.run(cmd + ["--", "python", "scripts/ensemble_components.py"], check=True, cwd=root)
