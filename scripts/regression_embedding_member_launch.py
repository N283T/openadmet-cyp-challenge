"""Launch the frozen-embedding member screen against the incumbent parent-cohort OOF."""

import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("embedding_run", help="shared run holding external-embeddings.parquet")
    parser.add_argument("--config", default="configs/regression/monroe-member-screen-v1.json")
    parser.add_argument(
        "--purpose", default="Screen frozen embedding readouts as regression members"
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    shared, regression = root / "artifacts/runs/shared", root / "artifacts/runs/regression"
    inputs = {
        "embeddings.parquet": shared / args.embedding_run / "outputs/external-embeddings.parquet",
        "chemeleon-embeddings.parquet": shared
        / "20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
        "test-mapping.parquet": shared
        / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
        "tabpfn.ckpt": shared
        / "20260910T122343854890Z-3a03f908/inputs/056-tabpfn-v2.6-regressor-v2.6_default.ckpt",
        "tabicl.ckpt": shared
        / "20260911T093102177307Z-beddd36f/outputs/tabicl-regressor-v2-20260212.ckpt",
        "control-1a2.parquet": regression
        / "20260918T034159381060Z-46314b50/outputs/1A2/add_dense_and_fast_pool/oof.parquet",
        "control-2c9.parquet": regression / "20260917T123601175420Z-e1560701/outputs/oof.parquet",
        "control-2d6.parquet": regression / "20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
        "control-3a4.parquet": regression / "20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    }
    command = ["pixi", "run", "cyp", "run", "--track", "regression", "--purpose", args.purpose]
    command += ["--config", str(root / args.config)]
    for alias, path in inputs.items():
        command += ["--input", f"{alias}={path}"]
    subprocess.run(
        command + ["--", "python", "scripts/regression_embedding_member.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
