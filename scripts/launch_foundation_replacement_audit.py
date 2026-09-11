"""Declare all updated model predictions and controls for final auditing."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = {
    "models": {},
    "database": {"host": str(root / ".runtime/socket"), "port": 5434, "dbname": "cyp_challenge"},
}
for backend in ["tabicl", "tabpfn"]:
    cfg["models"][backend] = []
    runs = re.findall(
        r"run_id=(\S+)", (root / f"artifacts/launches/replacement-{backend}.log").read_text()
    )
    assert len(runs) == 5
    for r in runs:
        out = root / "artifacts/runs/regression" / r / "outputs"
        suffix = "tabicl.joblib" if backend == "tabicl" else "tabpfn_fit"
        assert len(list(out.glob("*." + suffix))) == 4
        cfg["models"][backend].append(
            {
                "run": r,
                "summary": add(out / "summary.json"),
                "oof": add(out / "tabicl-parent-oof.parquet"),
                "test": add(out / "test-predictions.parquet"),
            }
        )
cfg["lgbm"] = add(
    root / "artifacts/runs/regression/20260911T130021181456Z-60011f3a/outputs/oof.parquet"
)
cfg["previous_tabicl"] = add(
    root / "artifacts/runs/regression/20260911T111952347762Z-e9f6c457/outputs/oof.parquet"
)
cfg["mapping"] = add(
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"
)
cfg["scorer"] = add(root / "artifacts/metric-reference/custom_scoring_functions.py")
p = root / "configs/parent-cv/foundation-replacement-audit.json"
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
    "Audit complete TabICL and TabPFN replacement including five-fold test predictions",
    "--config",
    str(p),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(
    cmd + ["--", "python", "scripts/audit_foundation_replacement.py"], check=True, cwd=root
)
