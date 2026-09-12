"""Launch the joint efficacy control and two public-loss weights."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
cfg = root / "configs/d6_joint_external.json"
with cfg.open("x") as f:
    json.dump(
        dict(
            public_weights=[0, 0.05, 0.30],
            max_epochs=80,
            patience=10,
            public_samples_per_epoch=4096,
            split_run="20260911T104703655337Z-a3923255",
        ),
        f,
        indent=2,
    )
inputs = [
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    b / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    b / "shared/20260912T101212878172Z-616de5d8/outputs/efficacy-targets.parquet",
    b / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    b / "shared/20260912T060426566774Z-7604ffeb/outputs/best.pt",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Joint D6 external efficacy graph training, matched no-public and publicweights .05/.30",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd.extend(["--input", str(p)])
subprocess.run(cmd + ["--", "python", "scripts/d6_joint_external.py"], cwd=root, check=True)
