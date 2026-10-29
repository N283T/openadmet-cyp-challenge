"""Launch the three-option D6 comparison in public-metric units."""

import argparse
import subprocess
from pathlib import Path

PROBE_RUN = "artifacts/runs/classification/20261003T042011147407Z-15430428/outputs"
PREREGISTRATION = "artifacts/tracked/cyp-d6-probe-preregistration-20261003T025354186349Z-621b6ad3"
INPUTS = {
    "blind_swaps.csv": f"{PROBE_RUN}/blind_swaps.csv",
    "marginal_precision.csv": f"{PREREGISTRATION}/marginal_precision.csv",
}
PURPOSE = (
    "Compare three D6 options on the recorded sizing basis of 374 scored rows and 26 scored "
    "positives: the gate-filtered vector at an unchanged 49 calls, the pre-registered monotone "
    "extension to the top 94 rows, and holding the deployed vector. Each option is simulated "
    "exactly rather than linearised, with a scored-fraction sensitivity and the probability "
    "that the extension yields a decisive outcome; read-only, no candidate or upload"
)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-option-comparison-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        PURPOSE,
        "--config",
        str(root / args.config),
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_option_comparison.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
