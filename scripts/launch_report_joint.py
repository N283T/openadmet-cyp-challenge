"""Record joint-feature preparation, matched LGBM comparison and independent audit."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs"
config_dir = root / "configs/report-joint"
config_dir.mkdir(exist_ok=True)


def launch(name, track, config, inputs, script, purpose):
    p = config_dir / (name + ".json")
    with p.open("x") as f:
        json.dump(config, f, indent=2)
    cmd = ["pixi", "run", "cyp", "run", "--track", track, "--purpose", purpose, "--config", str(p)]
    for f in inputs:
        cmd += ["--input", str(f)]
    r = subprocess.run(cmd + ["--", "python", script], cwd=root, text=True, capture_output=True)
    print(r.stdout, flush=True)
    if r.returncode:
        print(r.stderr)
    r.check_returncode()
    return re.search(r"run_id=(\S+)", r.stdout)[1]


features = launch(
    "features_v1",
    "shared",
    dict(feature_sets=["geometry", "substrate", "jazzy"]),
    [
        base / "shared/20260912T092823285236Z-3bc4bd9d/outputs/geometry-features.parquet",
        base / "shared/20260912T092919346028Z-bc6e15fc/outputs/substrate-features.parquet",
        base / "shared/20260912T085529851955Z-9df14c0c/outputs/jazzy-features.parquet",
    ],
    "scripts/report_joint_features.py",
    "Join geometry substrate and Jazzy features for matched combined control",
)
d = base / f"shared/{features}/outputs"
columns = json.loads((d / "feature-columns.json").read_text())
config = dict(
    split_run="20260911T104703655337Z-a3923255",
    backend="lgbm",
    control_tolerance=0,
    arms=dict(
        control=[],
        joint_geometry_substrate=columns["geometry"] + columns["substrate"],
        joint_geometry_substrate_jazzy=columns["geometry"]
        + columns["substrate"]
        + columns["jazzy"],
    ),
)
fit = launch(
    "paired_v1",
    "regression",
    config,
    [
        base / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
        base / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        base / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
        d / "joint-features.parquet",
        base / "regression/20260912T083010595992Z-1e7001de/outputs/external-oof.parquet",
    ],
    "scripts/report_geometry_readout.py",
    "Matched simultaneous geometry substrate and Jazzy additions",
)
launch(
    "audit_v1",
    "regression",
    dict(
        model_run=fit,
        database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
    ),
    [
        base / f"regression/{fit}/outputs/geometry-oof.parquet",
        base
        / "regression/20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
        base / "regression/20260912T072935321335Z-1e8d2643/outputs/external-oof.parquet",
        base / f"regression/{fit}/outputs/memberships.parquet",
        base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ],
    "scripts/report_geometry_audit.py",
    "Audit simultaneous geometry substrate and Jazzy OOF additions",
)
