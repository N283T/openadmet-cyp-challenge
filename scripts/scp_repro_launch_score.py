"""Score a recovered native member using the parent's fixed reference."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("recipe")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config = root / "configs/scp-repro" / f"score-{args.run}.json"
    config.write_text(json.dumps(dict(recipe=args.recipe, export_run=args.run), indent=2) + "\n")
    inputs = [
        root / "artifacts/metric-reference/custom_scoring_functions.py",
        root
        / "artifacts/runs/regression/20260912T094021204217Z-84c28266/outputs/global-oof.parquet",
        root / f"artifacts/runs/regression/{args.run}/outputs/scp-oof.parquet",
    ]
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Score completed native {args.recipe} member against fixed previous global OOF",
        "--config",
        str(config),
    ]
    for path in inputs:
        assert path.is_file()
        command += ["--input", str(path)]
    subprocess.run(
        command + ["--", "python", "scripts/scp_repro_score_available.py"], check=True, cwd=root
    )


if __name__ == "__main__":
    main()
