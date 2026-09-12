"""Run the two predeclared single-target interval adaptation controls sequentially."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--matching-logfc", action="store_true")
args = parser.parse_args()
tag = "matching_logfc" if args.matching_logfc else "single"
root = Path(__file__).resolve().parents[1]
data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
split = "20260911T104703655337Z-a3923255"
inputs = [
    root / "artifacts/runs/shared/20260912T060426566774Z-7604ffeb/outputs/best.pt",
    data / "compounds.parquet",
    data / "measurements.parquet",
    root / f"artifacts/runs/regression/{split}/outputs/parent-assignments.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
if args.matching_logfc:
    inputs.append(
        root / "artifacts/runs/shared/20260910T100337258849Z-927eac1c/outputs/targets.parquet"
    )
config_dir = root / "configs/focus-single-target"
config_dir.mkdir(exist_ok=True)
manifest = {}
for target, n in [("CYP1A2", 1372), ("CYP2D6", 1457)]:
    cfg = dict(
        endpoint=target + "_pIC50_direct_inhibition",
        split_run=split,
        expected_oof_rows=n,
        method="focus_" + tag + "_" + target,
        matching_logfc=args.matching_logfc,
        pIC50_weight=2.0,
        training_loss="pic50_interval",
        checkpoint_selection="pic50_st_rae",
        max_epochs=160,
        patience=20,
        num_threads=2,
        seed=42,
        inner_fraction=0.1,
    )
    path = config_dir / (target + "_" + tag + "_interval_v2.json")
    with path.open("x") as stream:
        json.dump(cfg, stream, indent=2)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Single-target interval adaptation control: " + target,
        "--config",
        str(path),
    ]
    for item in inputs:
        command += ["--input", str(item)]
    result = subprocess.run(
        command + ["--", "python", "scripts/focus_single_target.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    (root / f"artifacts/launches/focus_{tag}_{target}_v2.log").write_text(
        result.stdout + result.stderr
    )
    print(result.stdout[-2000:], flush=True)
    result.check_returncode()
    run = re.findall(r"run_id=(\S+)", result.stdout)
    assert len(run) == 1
    manifest[target] = run[0]
    (root / f"artifacts/launches/focus_{tag}_v2_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )
