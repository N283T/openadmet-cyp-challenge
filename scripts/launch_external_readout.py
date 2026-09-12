"""Launch one immutable readout plus matched audit for a completed feature asset."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--features", type=Path, required=True)
parser.add_argument("--method", required=True)
parser.add_argument("--backend", choices=["lgbm", "tabpfn"], required=True)
parser.add_argument("--prefix", default="emb_")
parser.add_argument("--pca", type=int, default=0)
parser.add_argument("--nested", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
assert re.fullmatch(r"[a-z0-9_]+", args.method)
files = [
    args.features.resolve(),
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet",
]
cfg = dict(
    method=args.method,
    backend=args.backend,
    feature_file=files[0].name,
    feature_prefix=args.prefix,
    pca=args.pca,
    split_run="20260911T104703655337Z-a3923255",
)
if args.nested:
    files.append(args.nested.resolve())
    cfg["nested_file"] = args.nested.name
if args.backend == "tabpfn":
    weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
    cfg["weights"] = f"{len(files):03d}-{weights.name}"
    files.append(weights)
config_dir = root / "configs/external-readout"
config_dir.mkdir(exist_ok=True)


def launch(name, config, inputs, script, purpose):
    p = config_dir / (name + ".json")
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
    for item in inputs:
        command.extend(["--input", str(item)])
    result = subprocess.run(
        command + ["--", "python", script], cwd=root, text=True, capture_output=True
    )
    (root / "artifacts/launches" / (name + ".log")).write_text(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()
    runs = re.findall(r"run_id=(\S+)", result.stdout)
    assert len(runs) == 1
    return runs[0]


r = launch(
    args.method,
    cfg,
    files,
    "scripts/external_embedding_readout.py",
    "Matched fixed-CV evaluation of external representation " + args.method,
)
baseline = (
    root / "artifacts/runs/regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet"
)
launch(
    args.method + "_audit",
    dict(
        model_run=r,
        database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
    ),
    [
        root / f"artifacts/runs/regression/{r}/outputs/external-oof.parquet",
        baseline,
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ],
    "scripts/audit_external_readout.py",
    "Matched scores and fixed equal blends for " + args.method,
)
