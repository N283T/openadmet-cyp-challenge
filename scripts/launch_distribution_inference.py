"""Launch each regression fold with remapped immutable inference inputs."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parent = root / "artifacts/runs/shared/20260910T132141471698Z-343eb5ec"
original = json.loads((parent / "config.json").read_text())
for fold in range(5):
    cfg = json.loads(json.dumps(original))
    cfg["folds"] = [cfg["folds"][fold]]
    cfg["folds"][0]["endpoints"] = {
        k: v for k, v in cfg["folds"][0]["endpoints"].items() if not k.endswith("is_TDI")
    }
    cfg["weights"].pop("classifier")
    cfg["quantiles"] = [0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9]
    cfg["scope"] = "Regression only; fixed distribution summaries, original model weights and fits"
    paths = []

    def remap(obj, paths=paths):
        if isinstance(obj, dict):
            return {k: remap(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [remap(v) for v in obj]
        if isinstance(obj, str) and (parent / "inputs" / obj).is_file():
            p = parent / "inputs" / obj
            if p not in paths:
                paths.append(p)
            return f"{paths.index(p):03d}-{p.name}"
        return obj

    cfg = remap(cfg)
    path = root / f"configs/calibration/distribution-fold{fold}-v1.json"
    with path.open("x") as stream:
        stream.write(json.dumps(cfg, indent=2) + "\n")
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"TabPFN 2.6 distribution readout fold {fold}; "
        "mean reload audit plus median/quantiles",
        "--config",
        str(path),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(
        cmd + ["--", "python", "scripts/distribution_inference.py"], check=True, cwd=root
    )
