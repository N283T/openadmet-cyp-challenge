"""Launch recorded geometry readout and fixed blend audit."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs"
features = base / "shared/20260912T092823285236Z-3bc4bd9d/outputs"
meta = json.loads((features / "feature-report.json").read_text())


def launch(tag, cfg, inputs, script, purpose):
    p = root / "configs/report_geometry" / f"{tag}.json"
    with p.open("x") as f:
        json.dump(cfg, f, indent=2)
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        purpose,
        "--config",
        str(p),
    ]
    for f in inputs:
        cmd += ["--input", str(f)]
    r = subprocess.run(cmd + ["--", "python", script], cwd=root, text=True, capture_output=True)
    (root / "artifacts/launches" / f"report_geometry_{tag}.log").write_text(r.stdout + r.stderr)
    print(r.stdout, r.stderr, flush=True)
    r.check_returncode()
    return re.search(r"run_id=(\S+)", r.stdout)[1]


fit = launch(
    "pfn_paired",
    dict(
        split_run="20260911T104703655337Z-a3923255",
        backend="tabpfn",
        control_tolerance=1e-5,
        arms={
            "control": [],
            "geometry_mean_std": meta["topology_columns"]
            + meta["mean_columns"]
            + meta["std_columns"],
        },
    ),
    [
        base / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
        base / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        base / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
        features / "geometry-features.parquet",
        base / "regression/20260912T083032513346Z-568bc526/outputs/external-oof.parquet",
        Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt",
    ],
    "scripts/report_geometry_readout.py",
    "Controlled ligand-only geometry additions to frozen600 plus nested log2fc",
)
launch(
    "pfn_audit",
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
    "Audit geometry control and fixed 25/50 percent ensemble additions",
)
