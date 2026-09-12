"""Declare native arrays and recovered exports for independent recorded audit."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("training_run")
    parser.add_argument("export_run")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs"
    raw = runs / "regression" / args.training_run / "outputs"
    export = runs / "regression" / args.export_run / "outputs"
    original = json.loads((raw / "effective-config.json").read_text())["config"]
    base = runs / "shared/20260909T132734119233Z-a7063ccd/outputs"
    paths = dict(
        labels=base / "measurements.parquet",
        parent=runs
        / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        assignment=raw / "assignments.parquet",
        mapping=base / "test_mapping.parquet",
        oof=export / "scp-oof.parquet",
        test=export / "test-predictions.parquet",
        fold_test=export / "fold-test-predictions.parquet",
    )
    for fold in range(5):
        paths[f"val{fold}"] = raw / f"fold{fold}/all-head-validation.parquet"
        paths[f"test{fold}"] = raw / f"fold{fold}/all-head-test.parquet"
        paths[f"member{fold}"] = raw / f"fold{fold}/membership.parquet"
    cfg = {key: f"{i:03d}-{path.name}" for i, (key, path) in enumerate(paths.items())}
    cfg.update(
        training_run=args.training_run,
        export_run=args.export_run,
        multi=original["recipe"] in ["union18", "aux8"],
        source_cv=original["cv_policy"] == "source_scaffold",
    )
    config = root / "configs/scp-repro" / f"audit-{args.export_run}.json"
    config.write_text(json.dumps(cfg, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Independent canonical membership and exact native-array export audit",
        "--config",
        str(config),
    ]
    for path in paths.values():
        assert path.is_file(), path
        command += ["--input", str(path)]
    subprocess.run(
        command + ["--", "python", "scripts/scp_repro_export_audit.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
