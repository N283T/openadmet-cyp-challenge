"""Declare source matrices and all saved boosters for independent classification audit."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification/20260912T115220209511Z-76facc71"
inputs = [run / "config.json", *sorted((run / "inputs").iterdir())]
outputs = ["oof.parquet", "fold-test-probabilities.parquet", "membership.parquet", "models.json"]
inputs += [run / "outputs" / name for name in outputs]
inputs += sorted((run / "outputs").glob("*.txt"))
assert len(list((run / "outputs").glob("*.txt"))) == 240
assert (run / "result.json").is_file()
inputs += [run / "result.json"]
config = root / "configs/classification/lgbm-independent-audit.json"
config.write_text(
    json.dumps(
        dict(
            source_config="000-config.json",
            source_run=run.name,
            mode="Read-only reconstruction of all120innerand120outerboosters; no fitting",
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
    "classification",
    "--purpose",
    "Independently recover inner predictions and audit all LGBM OOF/test probabilities "
    "and heldout MCC thresholds",
    "--config",
    str(config),
]
for p in inputs:
    assert p.is_file(), p
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/classification_lgbm_audit.py"], check=True, cwd=root)
