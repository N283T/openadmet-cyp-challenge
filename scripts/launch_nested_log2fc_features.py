"""Launch auditable feature assembly from the adopted ten auxiliary pair models."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
reg_path = root / "configs/regression/current-log2fc.json"
readout_path = root / "configs/regression/current-parent-readout.json"
reg = json.loads(reg_path.read_text())
readout = json.loads(readout_path.read_text())
paths = []


def add(path):
    paths.append(path)
    return f"{len(paths) - 1:03d}-{path.name}"


cfg = {
    "registry": add(reg_path),
    "readout_registry": add(readout_path),
    "parent_assignments": add(
        root
        / "artifacts/runs/regression"
        / readout["split_run"]
        / "outputs/parent-assignments.parquet"
    ),
    "aux_assignments": add(
        root / "artifacts/runs/regression" / reg["split_run"] / "outputs/assignments.parquet"
    ),
    "models": [],
    "dtype": "float32",
    "population": "6146 challenge training compounds",
    "routing": (
        "Recipient fold r != outer o uses model excluding o,r; "
        "r=o mean of four models excluding o."
    ),
}
for model in reg["models"]:
    base = root / "artifacts/runs/shared" / model["run"] / "outputs"
    cfg["models"].append(
        {
            **model,
            "summary": add(base / "summary.json"),
            "membership": add(base / "membership.parquet"),
            "predictions": [add(p) for p in sorted(base.glob("predictions-*.parquet"))],
        }
    )
config_path = root / "configs/external-readout/nested-log2fc.json"
config_path.parent.mkdir(parents=True, exist_ok=True)
if config_path.exists():
    assert json.loads(config_path.read_text()) == cfg
else:
    with config_path.open("x") as f:
        json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Assemble nested log2fc features from ten pair models with full membership audit",
    "--config",
    str(config_path),
]
for path in paths:
    cmd.extend(["--input", str(path)])
subprocess.run(
    cmd + ["--", "python", "scripts/prepare_nested_log2fc_features.py"], check=True, cwd=root
)
