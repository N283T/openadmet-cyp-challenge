"""Launch fold-wise MONROE fine-tuning against the incumbent parent-cohort OOF."""

import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph_run", help="shared run holding graphs-seed*.pkl")
    parser.add_argument("--config", required=True)
    parser.add_argument("--purpose", required=True)
    parser.add_argument("--init-config", help="checkpoint config.json override")
    parser.add_argument("--init-weights", help="checkpoint weights.pt override")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    shared, regression = root / "artifacts/runs/shared", root / "artifacts/runs/regression"
    checkpoint = Path.home() / "ghq/github.com/blazejba/monroe/checkpoint"
    inputs = {
        f"graphs-seed{s}.pkl": shared / args.graph_run / f"outputs/graphs-seed{s}.pkl"
        for s in range(1, 6)
    }
    inputs.update(
        {
            "monroe-config.json": Path(args.init_config or checkpoint / "config.json"),
            "monroe-weights.pt": Path(args.init_weights or checkpoint / "weights.pt"),
            "measurements.parquet": shared
            / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
            "aliases.parquet": shared / "20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
            "primary.csv": root / "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
            "parent-assignments.parquet": regression
            / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
            "test-mapping.parquet": shared
            / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
            "custom_scoring_functions.py": regression
            / "20260917T121017505228Z-433d9acb/inputs/custom_scoring_functions.py",
            "tabpfn.ckpt": shared
            / "20260910T122343854890Z-3a03f908/inputs/056-tabpfn-v2.6-regressor-v2.6_default.ckpt",
            "control-1a2.parquet": regression
            / "20260918T034159381060Z-46314b50/outputs/1A2/add_dense_and_fast_pool/oof.parquet",
            "control-2c9.parquet": regression
            / "20260917T123601175420Z-e1560701/outputs/oof.parquet",
            "control-2d6.parquet": regression
            / "20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
            "control-3a4.parquet": regression
            / "20260917T130718808881Z-dfa02730/outputs/oof.parquet",
        }
    )
    command = ["pixi", "run", "cyp", "run", "--track", "regression", "--purpose", args.purpose]
    command += ["--config", str(root / args.config)]
    for alias, path in inputs.items():
        command += ["--input", f"{alias}={path}"]
    python = str(root / ".pixi/envs/monroe/bin/python")
    subprocess.run(command + ["--", python, "scripts/monroe_finetune.py"], cwd=root, check=True)


if __name__ == "__main__":
    main()
