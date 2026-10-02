"""Record the independent restart comparison without historical-model inputs."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-zero-start-report-v1.json")
    args = parser.parse_args()
    config = root / args.config
    cfg = json.loads(config.read_text())
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Compare independent D6 supervision and cross-enzyme learning "
        "with paired scaffold bootstrap",
        "--config",
        str(config),
    ]
    for version, run in cfg["sources"].items():
        for name in ["summary.json", "oof.parquet"]:
            path = root / "artifacts/runs/classification" / run / "outputs" / name
            if not path.is_file():
                raise FileNotFoundError(path)
            command.extend(["--input", f"{version}-{name}={path}"])
    assignments = (
        root
        / "artifacts/runs/classification"
        / cfg["sources"]["v1"]
        / "outputs/assignments.parquet"
    )
    command.extend(["--input", f"assignments.parquet={assignments}"])
    for alias, run in cfg.get("extra_assignments", {}).items():
        path = root / "artifacts/runs/classification" / run / "outputs/assignments.parquet"
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/report_cyp_d6_zero_start.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
