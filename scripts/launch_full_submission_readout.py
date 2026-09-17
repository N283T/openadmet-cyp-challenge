"""Launch recorded full-data readout fits, retaining all immutable inputs."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--method", required=True)
p.add_argument("--features", type=Path, required=True)
p.add_argument("--prefix", default="emb_")
p.add_argument("--extra-prefix", action="append", default=[])
p.add_argument("--dimensions", type=int, required=True)
p.add_argument("--backend", choices=["tabpfn", "tabicl", "tuned_lgbm", "fixed_lgbm"], required=True)
p.add_argument("--endpoint", action="append", default=[])
p.add_argument("--auxiliary", type=Path)
p.add_argument("--pca", type=int, default=0)
p.add_argument("--tabpfn-version", default="v2.6")
p.add_argument("--tabpfn-weights", type=Path)
args = p.parse_args()
assert re.fullmatch(r"[a-z0-9_]+", args.method)
inputs = {}


def add(alias, path):
    path = Path(path).resolve()
    assert path.is_file(), path
    assert alias not in inputs
    inputs[alias] = path
    return alias


data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = dict(
    method=args.method,
    backend=args.backend,
    prefix=args.prefix,
    prefixes=[args.prefix, *args.extra_prefix],
    expected_dimensions=args.dimensions,
    pca=args.pca,
    features=add("features.parquet", args.features),
    measurements=add("measurements.parquet", data / "measurements.parquet"),
    assignments=add(
        "assignments.parquet",
        root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet",
    ),
    aliases=add("aliases.parquet", data / "aliases.parquet"),
    reference=add("reference.csv", root / "data/raw/cyp-challenge-TEST-BLINDED.csv"),
)
if args.endpoint:
    cfg["endpoints"] = args.endpoint
if args.auxiliary:
    cfg["auxiliary"] = add("auxiliary.parquet", args.auxiliary)
if args.backend == "tabpfn":
    weights = args.tabpfn_weights
    if weights is None:
        weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
    cfg["weights"] = add("weights.safetensors", weights)
    cfg["tabpfn_version"] = args.tabpfn_version
    cfg["tabpfn_estimators"] = 8
elif args.backend == "tabicl":
    cfg["weights"] = add(
        "weights.ckpt",
        root
        / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
        / "tabicl-regressor-v2-20260212.ckpt",
    )
elif args.backend == "tuned_lgbm":
    cfg["selection_summary"] = add(
        "selection-summary.json",
        root / "artifacts/runs/regression/20260912T072056445264Z-e6b6cbcb/outputs/summary.json",
    )
    cfg["deployment_selection"] = (
        "Refit each of five inner-CV-selected parameter sets on all endpoint labels, "
        "equal mean within component."
    )
directory = root / "configs/submission-full"
directory.mkdir(exist_ok=True)
config = directory / f"{args.method}.json"
if config.exists():
    assert json.loads(config.read_text()) == cfg, config
else:
    with config.open("x") as f:
        json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Full-data refit for regression candidate component: " + args.method,
    "--config",
    str(config),
]
for alias, path in inputs.items():
    cmd += ["--input", f"{alias}={path}"]
result = subprocess.run(
    cmd + ["--", "python", "scripts/full_submission_readout.py"],
    cwd=root,
    text=True,
    capture_output=True,
)
(root / "artifacts/launches" / f"{args.method}.log").write_text(result.stdout + result.stderr)
print(result.stdout, flush=True)
print(result.stderr, flush=True)
result.check_returncode()
