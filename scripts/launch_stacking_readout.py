"""Launch recorded PCA + nested log2fc readouts from completed auxiliary runs."""

import argparse
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/shared"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fold", type=int, action="append")
    args = parser.parse_args()
    auxiliary = {}
    for run in re.findall(
        r"run_id=(\S+)", (ROOT / "artifacts/launches/stacking-aux-20260910.log").read_text()
    ):
        out = RUNS / run / "outputs"
        if not (out / "summary.json").exists():
            continue
        summary = json.loads((out / "summary.json").read_text())
        auxiliary[tuple(summary["config"]["excluded_folds"])] = (run, out)
    for fold in args.fold if args.fold is not None else range(5):
        assert sum(fold in pair for pair in auxiliary) == 4
        config = json.loads(
            (ROOT / f"configs/readout/tabpfn-embedding-fold{fold}.json").read_text()
        )
        config.update(
            features="foundation_pca256_plus_nested_log2fc4",
            pca_components=256,
            foundation_run="20260910T012630927160Z-e9576045",
            auxiliary=[],
        )
        files = []

        def add(path, files=files):
            index = len(files)
            files.append(path)
            return f"{index:03d}-{path.name}"

        for run, name in [
            (config["dataset_run"], "measurements.parquet"),
            (config["dataset_run"], "aliases.parquet"),
            (config["split_run"], "assignments.parquet"),
        ]:
            add(RUNS / run / "outputs" / name)
        for kind in ["regressor", "classifier"]:
            add(Path.home() / f".cache/tabpfn/tabpfn-v2.6-{kind}-v2.6_default.ckpt")
        for path in sorted(
            (RUNS / config["foundation_run"] / "outputs").glob("chemeleon-*.parquet")
        ):
            add(path)
        for pair, (run, out) in sorted(auxiliary.items()):
            if fold not in pair:
                continue
            config["auxiliary"].append(
                dict(
                    run=run,
                    excluded_folds=list(pair),
                    summary=add(out / "summary.json"),
                    membership=add(out / "membership.parquet"),
                    predictions=[add(p) for p in sorted(out.glob("predictions-*.parquet"))],
                )
            )
        config_path = ROOT / f"configs/stacking/readout-fold{fold}.json"
        assert not config_path.exists(), "Do not overwrite an existing campaign config"
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        command = [
            "pixi",
            "run",
            "cyp",
            "run",
            "--track",
            "shared",
            "--purpose",
            f"UMAP fold {fold}: foundation PCA256 plus nested log2fc4 to TabPFN-2.6",
            "--config",
            str(config_path),
        ]
        for p in files:
            command += ["--input", str(p)]
        subprocess.run(
            command + ["--", "python", "scripts/stacking_readout.py"], check=True, cwd=ROOT
        )


if __name__ == "__main__":
    main()
