"""Record an independent audit of a completed focused multitask run."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--logfc-weight", type=float, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = root / "artifacts/runs/regression" / args.run / "outputs"
    config = root / "configs/external-finetune" / f"{args.run}_focus_verify.json"
    with config.open("x") as stream:
        json.dump(dict(model_run=args.run, logfc_weight=args.logfc_weight), stream, indent=2)
    inputs = [output / "summary.json", output / "finetune-oof.parquet"]
    for fold in range(5):
        inputs += [output / f"fold{fold}/membership.parquet", output / f"fold{fold}/history.jsonl"]
    assert all(path.is_file() for path in inputs)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Independent focused multitask membership, checkpoint selection and OOF verification",
        "--config",
        str(config),
    ]
    for path in inputs:
        command.extend(["--input", str(path)])
    subprocess.run(
        command + ["--", "python", "scripts/verify_focus_multitask.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
