"""Run a fixed objective comparison and downstream fits, with immutable run records."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--dataset-run", required=True)
parser.add_argument("--tag", required=True)
parser.add_argument(
    "--modes", nargs="+", default=["classification", "ranking", "continuous", "mixed"]
)
parser.add_argument("--epochs", type=int, default=60)
parser.add_argument("--steps", type=int, default=200)
parser.add_argument("--conditioned", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
assert re.fullmatch(r"[a-z0-9_]+", args.tag)
folder = root / f"artifacts/runs/shared/{args.dataset_run}/outputs"
files = [
    folder / "external-compounds.parquet",
    folder / "observations.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
]
assert all(p.is_file() for p in files)
configs = root / "configs/external-training"
configs.mkdir(exist_ok=True)
for mode in args.modes:
    assert mode in ["continuous", "ranking", "classification", "mixed"]
    name = args.tag + "_" + mode
    cfg = dict(
        mode=mode,
        dataset_run=args.dataset_run,
        seed=42,
        max_epochs=args.epochs,
        patience=10,
        batch_size=128,
        steps_per_epoch=args.steps,
        min_train_observations=5,
        rank_min_gap=0.25,
        initialization="random",
        width=600,
        depth=5,
        sampling="uniform objective then source then assay",
        purpose="Controlled external objective comparison on shared dataset and seed",
    )
    if args.conditioned:
        cfg.update(conditioning_group="cyp_measurement", conditioning_variant="affine")
    path = configs / (name + ".json")
    with path.open("x") as stream:
        json.dump(cfg, stream, indent=2)
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "External600 pretraining with " + mode + " objective; " + args.tag,
        "--config",
        str(path),
    ]
    for p in files:
        cmd += ["--input", str(p)]
    result = subprocess.run(
        cmd
        + [
            "--",
            "python",
            "scripts/train_external_conditioned.py"
            if args.conditioned
            else "scripts/train_external_assay.py",
        ],
        cwd=root,
        capture_output=True,
        text=True,
    )
    (root / "artifacts/launches" / (name + ".log")).write_text(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()
    runs = re.findall(r"run_id=(\S+)", result.stdout)
    assert len(runs) == 1
    feature = root / f"artifacts/runs/shared/{runs[0]}/outputs/external-embeddings.parquet"
    for backend in ["lgbm", "tabpfn"]:
        subprocess.run(
            [
                "pixi",
                "run",
                "python",
                "scripts/launch_external_readout.py",
                "--features",
                str(feature),
                "--method",
                name + "_" + backend,
                "--backend",
                backend,
            ],
            cwd=root,
            check=True,
        )
