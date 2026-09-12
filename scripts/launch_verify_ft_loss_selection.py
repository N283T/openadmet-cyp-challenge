"""Record loss-selection smoke verification assets."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/regression/20260912T060513219458Z-97152a22/outputs"
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = dict(scorer=add(root / "artifacts/metric-reference/custom_scoring_functions.py"), folds=[])
for fold in range(5):
    b = base / f"fold{fold}"
    cfg["folds"].append(
        dict(
            fold=fold,
            history=add(b / "history.jsonl"),
            checkpoint=add(b / "best.pt"),
            membership=add(b / "membership.parquet"),
        )
    )
config = root / "configs/external-finetune/mae-st-selection-verify.json"
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
    "Verify inner ST checkpoint decisions and default/MAE loss numerics",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/verify_ft_loss_selection.py"], check=True, cwd=root)
