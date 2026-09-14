"""Launch the fixed-panel TDI structural-alert augmentation experiment."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    candidate = root / path
    assert candidate.is_file(), candidate
    paths.append(candidate)
    return f"{len(paths) - 1:03d}-{candidate.name}"


base = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = {
    name: add(f"{base}/{filename}.parquet")
    for name, filename in [
        ("compounds", "compounds"),
        ("labels", "measurements"),
        ("test_mapping", "test_mapping"),
    ]
}
cfg["assignments"] = add(
    "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
)
cfg["split_run"] = "20260909T135130268493Z-edba6fc0"
cfg["mordred"] = [
    add(path.relative_to(root))
    for path in sorted(
        (root / "artifacts/runs/shared/20260909T132918829974Z-bc81af58/outputs").glob(
            "mordred-*.parquet"
        )
    )
]
cfg["alert_reference"] = add(
    "references/official/cyp/structural-alerts/20260907T135802Z/snapshot.md"
)
cfg["endpoints"] = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
config = root / "configs/classification-round/tdi-alerts-v1.json"
config.write_text(json.dumps(cfg, indent=2) + "\n")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Test a fixed 12-category TDI structural-alert panel added to "
        "Morgan and Morgan+Mordred LightGBM"
    ),
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_tdi_alerts.py"],
    cwd=root,
    check=True,
)
