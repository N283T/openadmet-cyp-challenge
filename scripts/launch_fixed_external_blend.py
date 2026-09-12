"""Launch a prespecified equal blend and its official fixed-CV audit."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--tag", required=True)
parser.add_argument(
    "--component", nargs=2, action="append", required=True, metavar=("PARQUET", "COLUMN")
)
args = parser.parse_args()
assert re.fullmatch(r"[a-z0-9_]+", args.tag)
root = Path(__file__).resolve().parents[1]
folder = root / "configs/external-ensemble"
inputs = [Path(p).resolve() for p, _ in args.component]
config = dict(
    method=args.tag,
    split_run="20260911T104703655337Z-a3923255",
    components=[
        dict(file=f"{i:03d}-{path.name}", column=column)
        for i, (path, (_, column)) in enumerate(zip(inputs, args.component, strict=True))
    ],
    weight_policy="Fixed equal mean, no OOF fitting",
)


def launch(name, cfg, files, script):
    config_path = folder / (name + ".json")
    with config_path.open("x") as stream:
        json.dump(cfg, stream, indent=2)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Fixed equal blend or matched official audit: " + name,
        "--config",
        str(config_path),
    ]
    for path in files:
        assert path.is_file(), path
        command += ["--input", str(path)]
    result = subprocess.run(
        command + ["--", "python", script], cwd=root, capture_output=True, text=True
    )
    with (root / "artifacts/launches" / (name + ".log")).open("x") as stream:
        stream.write(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()
    return re.findall(r"run_id=(\S+)", result.stdout)[0]


run = launch(args.tag, config, inputs, "scripts/fixed_external_blend.py")
launch(
    args.tag + "_audit",
    dict(
        model_run=run,
        database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
    ),
    [
        root / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet",
        root
        / "artifacts/runs/regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ],
    "scripts/audit_external_readout.py",
)
