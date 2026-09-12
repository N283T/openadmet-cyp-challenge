"""Launch paired Jazzy readout and audit with all immutable assets declared."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs"
parser = argparse.ArgumentParser()
parser.add_argument("--features", required=True, type=Path)
parser.add_argument("--mean-columns", nargs="+", required=True)
parser.add_argument("--std-columns", nargs="+", required=True)
parser.add_argument("--tag", default="jazzy_dg_v1")
parser.add_argument("--backend", choices=["lgbm", "tabpfn"], default="lgbm")
args = parser.parse_args()
assert re.fullmatch("[a-z0-9_]+", args.tag)
config_dir = root / "configs/jazzy"
config_dir.mkdir(exist_ok=True)


def launch(suffix, config, inputs, script, purpose):
    p = config_dir / f"{args.tag}_{suffix}.json"
    with p.open("x") as f:
        json.dump(config, f, indent=2)
    command = [
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
        command.extend(["--input", str(f)])
    r = subprocess.run(command + ["--", "python", script], cwd=root, text=True, capture_output=True)
    (root / "artifacts/launches" / f"{args.tag}_{suffix}.log").write_text(r.stdout + r.stderr)
    print(r.stdout, r.stderr, flush=True)
    r.check_returncode()
    return re.search(r"run_id=(\S+)", r.stdout)[1]


control_run = (
    "20260912T083010595992Z-1e7001de"
    if args.backend == "lgbm"
    else "20260912T083032513346Z-568bc526"
)
inputs = [
    base / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    base / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    base / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    args.features.resolve(),
    base / f"regression/{control_run}/outputs/external-oof.parquet",
]
if args.backend == "tabpfn":
    inputs.append(Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt")
fit = launch(
    "paired",
    dict(
        split_run="20260911T104703655337Z-a3923255",
        backend=args.backend,
        control_tolerance=0 if args.backend == "lgbm" else 1e-5,
        arms={
            "control": [],
            "jazzy_mean": args.mean_columns,
            "jazzy_mean_std": args.mean_columns + args.std_columns,
        },
    ),
    inputs,
    "scripts/jazzy_readout.py",
    "Controlled DG Jazzy additions to frozen600 and new nested log2fc",
)
launch(
    "audit",
    dict(
        model_run=fit,
        database=dict(
            host=str(root / ".runtime/socket"),
            port=5434,
            dbname="cyp_challenge",
        ),
    ),
    [
        base / f"regression/{fit}/outputs/jazzy-oof.parquet",
        base
        / "regression/20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
        base / "regression/20260912T072935321335Z-1e8d2643/outputs/external-oof.parquet",
        base / f"regression/{fit}/outputs/memberships.parquet",
        base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ],
    "scripts/audit_jazzy_readout.py",
    "Verify Jazzy paired controls and fixed25/50 ensemble additions",
)
