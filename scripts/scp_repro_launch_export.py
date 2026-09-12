"""Record a lossless canonical prediction export for a native training run."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    run = root / "artifacts/runs/regression" / args.run
    output = run / "outputs"
    original = json.loads((output / "effective-config.json").read_text())["config"]
    base = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    paths = dict(
        labels=base / "measurements.parquet",
        assignments=output / "assignments.parquet",
        parent=root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet",
        prediction_compounds=root
        / "artifacts/runs/shared/20260912T104336148000Z-ae7e37fc/outputs"
        / "prediction-compounds.parquet",
        test_mapping=base / "test_mapping.parquet",
    )
    for fold in range(5):
        paths[f"val{fold}"] = output / f"fold{fold}/all-head-validation.parquet"
        paths[f"test{fold}"] = output / f"fold{fold}/all-head-test.parquet"
    cfg = {key: f"{i:03d}-{path.name}" for i, (key, path) in enumerate(paths.items())}
    cfg.update(
        source_run=args.run,
        recipe=original["recipe"],
        cv_policy=original["cv_policy"],
        split_id=original["split_run"],
    )
    config = root / "configs/scp-repro" / f"export-{args.run}.json"
    config.write_text(json.dumps(cfg, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Recover canonical keyed OOF and five-fold test predictions "
        "from immutable native Chemprop outputs",
        "--config",
        str(config),
    ]
    for path in [*paths.values(), output / "summary.json"]:
        assert path.is_file(), path
        command += ["--input", str(path)]
    subprocess.run(command + ["--", "python", "scripts/scp_repro_export.py"], check=True, cwd=root)


if __name__ == "__main__":
    main()
