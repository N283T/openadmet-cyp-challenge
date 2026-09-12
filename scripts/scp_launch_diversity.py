"""Record the completed evaluator's fixed prediction diversity diagnostic."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("evaluation_run")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = (
        root / "artifacts/runs/regression" / args.evaluation_run / "outputs/all-model-oof.parquet"
    )
    assert path.is_file()
    cfg = root / "configs/scp-repro" / f"diversity-{args.evaluation_run}.json"
    cfg.write_text(
        json.dumps(
            dict(oof="000-all-model-oof.parquet", evaluation_run=args.evaluation_run), indent=2
        )
    )
    subprocess.run(
        [
            "pixi",
            "run",
            "cyp",
            "run",
            "--track",
            "regression",
            "--purpose",
            "Diagnose fixed D6 prediction/residual diversity and bias without fitting weights",
            "--config",
            str(cfg),
            "--input",
            str(path),
            "--",
            "python",
            "scripts/scp_diversity.py",
        ],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
