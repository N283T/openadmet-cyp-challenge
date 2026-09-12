"""Launch native published Chemprop recipes with immutable source/data inputs."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("recipe", choices=["aux8", "union18", "isoform6", "single1"])
    parser.add_argument("--source-cv", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = root / "artifacts/runs/shared/20260912T104004906145Z-97b9c270/outputs"
    data = root / "artifacts/runs/shared/20260912T104336148000Z-ae7e37fc/outputs"
    base = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    commit = "4cf07a8cb4366eb29111fb42396336b20b69ada0"
    split = "20260911T104703655337Z-a3923255"
    key = dict(aux8="challenge8", union18="union18", isoform6="d6_isoform6", single1="d6_single1")[
        args.recipe
    ]
    contract = json.loads((data / "contract.json").read_text())["models"][key]
    declared = dict(
        source_helpers=source / "src/workbench/endpoints/chemprop_utils.py",
        source_core=source / "src/workbench/training/chemprop_core.py",
        source_splits=root
        / "artifacts/runs/shared/20260912T104213710169Z-9388c902/outputs"
        / commit
        / "src/workbench/training/splits.py",
        source_template=source / "chemprop.template",
        table_file=data / contract["file"],
        assignments_file=root
        / f"artifacts/runs/regression/{split}/outputs/parent-assignments.parquet",
        connectivity_folds_file=data / "challenge-connectivity-folds.parquet",
        labels_file=base / "measurements.parquet",
        test_file=data / "prediction-compounds.parquet",
        test_mapping_file=base / "test_mapping.parquet",
    )
    config = {key: f"{i:03d}-{path.name}" for i, (key, path) in enumerate(declared.items())}
    config.update(
        recipe=args.recipe,
        heads=contract["heads"],
        n_primary=4 if args.recipe in ["aux8", "union18"] else 1,
        scored_endpoints=["CYP3A4", "CYP2C9", "CYP2D6", "CYP1A2"]
        if args.recipe in ["aux8", "union18"]
        else ["CYP2D6"],
        split_run="scp_source_scaffold_42" if args.source_cv else split,
        source_commit=commit,
        num_threads=2,
        num_workers=0,
        cv_policy="source_scaffold" if args.source_cv else "parent_adapted",
    )
    config["scored_endpoints"] = [
        target + "_pIC50_direct_inhibition" for target in config["scored_endpoints"]
    ]
    folder = root / "configs/scp-repro"
    folder.mkdir(exist_ok=True)
    path = folder / (args.recipe + ("-sourcecv" if args.source_cv else "") + ".json")
    path.write_text(json.dumps(config, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Published Workbench native Chemprop reproduction: {key}; {config['cv_policy']}; "
        "actual source defaults and five-fold prediction ensemble",
        "--config",
        str(path),
    ]
    for item in [*declared.values(), data / "contract.json"]:
        assert item.is_file(), item
        command += ["--input", str(item)]
    subprocess.run(command + ["--", "python", "scripts/scp_repro_train.py"], cwd=root, check=True)


if __name__ == "__main__":
    main()
