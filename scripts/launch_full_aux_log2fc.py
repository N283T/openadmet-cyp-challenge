"""Launch the full auxiliary refit with immutable inputs."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
inputs = []


def add(path):
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


registry = root / "configs/regression/current-log2fc.json"
reg = json.loads(registry.read_text())
cfg = {
    "registry": add(registry),
    "compounds": add(shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"),
    "targets": add(shared / "20260910T100337258849Z-927eac1c/outputs/targets.parquet"),
    "nested": add(shared / "20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"),
    "checkpoint": add(shared / "20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt"),
    "provenance": add(shared / "20260910T012612848389Z-bfeed434/outputs/provenance.json"),
    "pair_summaries": [add(shared / m["run"] / "outputs/summary.json") for m in reg["models"]],
    "nested_membership_audit": add(
        shared / "20260912T052549384296Z-bac5fb3a/outputs/membership-audit.json"
    ),
}
path = root / "configs/full-submission/full-aux-log2fc-v1.json"
path.parent.mkdir(parents=True, exist_ok=True)
if path.exists():
    assert json.loads(path.read_text()) == cfg
else:
    path.write_text(json.dumps(cfg, indent=2) + "\n")
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Full-data log2fc refit and cross-fitted train auxiliaries for submission",
    "--config",
    str(path),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/full_aux_log2fc.py"], cwd=root, check=True)
