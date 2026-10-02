"""Record complete historical ensemble refits under the new nested protocol."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-legacy-whole-nested-v1.json")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs"
    cfg_path = root / args.config
    cfg = json.loads(cfg_path.read_text())
    shared = runs / "shared/20260909T132734119233Z-a7063ccd/outputs"
    original = runs / "classification/20260912T115220209511Z-76facc71"
    inputs = {
        name: shared / name
        for name in ["compounds.parquet", "measurements.parquet", "aliases.parquet"]
    }
    inputs["assignments.parquet"] = (
        runs / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
    )
    inputs["class-embeddings.parquet"] = original / "inputs/004-external-embeddings.parquet"
    inputs["mixed-embeddings.parquet"] = original / "inputs/005-external-embeddings.parquet"
    inputs["external-best.pt"] = (
        runs / "classification/20260912T120135164747Z-1da37202/inputs/005-best.pt"
    )
    inputs["tabpfn-v2.6.ckpt"] = (
        runs
        / "classification/20260912T115242271820Z-51023dac/inputs"
        / "003-004-tabpfn-v2.6-classifier-v2.6_default.ckpt"
    )
    for i in range(cfg["mordred_files"]):
        paths = list((original / "inputs").glob(f"*-mordred-{i:05d}.parquet"))
        assert len(paths) == 1
        inputs[f"mordred-{i:05d}.parquet"] = paths[0]
    for path in [cfg_path, *inputs.values()]:
        if not path.is_file():
            raise FileNotFoundError(path)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        cfg["scientific_purpose"],
        "--config",
        str(cfg_path),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_zero_start.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
