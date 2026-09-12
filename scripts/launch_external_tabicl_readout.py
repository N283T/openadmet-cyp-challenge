"""Launch recorded TabICL readouts of frozen or five fold-specific embeddings."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--features", type=Path)
    source.add_argument("--fold-features", nargs=5, type=Path)
    source.add_argument("--finetune-dir", type=Path)
    parser.add_argument("--method", required=True)
    parser.add_argument("--nested", type=Path)
    parser.add_argument("--pca", type=int, default=256)
    parser.add_argument("--prefix", default="emb_")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    assert re.fullmatch(r"[a-z0-9_]+", args.method)
    assert 1 <= args.threads <= 2 and args.pca > 0
    root = Path(__file__).resolve().parents[1]
    split_run = "20260911T104703655337Z-a3923255"
    data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    assignments = root / f"artifacts/runs/regression/{split_run}/outputs/parent-assignments.parquet"
    weights = (
        root
        / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
        / "tabicl-regressor-v2-20260212.ckpt"
    )
    files = []

    def register(path):
        path = path.resolve()
        assert path.is_file(), path
        if path not in files:
            files.append(path)
        return f"{files.index(path):03d}-{path.name}"

    config = dict(
        method=args.method,
        split_run=split_run,
        assignments=register(assignments),
        measurements=register(data / "measurements.parquet"),
        weights=register(weights),
        fold_specific=args.features is None,
        prepare_only=args.prepare_only,
        pca=args.pca,
        threads=args.threads,
        feature_prefix=args.prefix,
        n_estimators=8,
        device="cuda",
        expected_oof_rows=6384,
    )
    if args.features:
        features = [args.features] * 5
    elif args.finetune_dir:
        features = [
            args.finetune_dir / f"fold{fold}/external-embeddings.parquet" for fold in range(5)
        ]
    else:
        features = args.fold_features
    config["feature_files"] = {str(fold): register(path) for fold, path in enumerate(features)}
    if args.nested:
        config["nested_file"] = register(args.nested)
    config_dir, log_dir = root / "configs/external-tabicl", root / "artifacts/launches"
    config_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)
    config_path = config_dir / f"{args.method}.json"
    log_path = log_dir / f"{args.method}.log"
    assert not config_path.exists() and not log_path.exists()
    with config_path.open("x") as stream:
        json.dump(config, stream, indent=2)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Fixed parent CV TabICL: embedding-only training PCA, then nested log2fc; "
        + ("CPU preparation only: " if args.prepare_only else "fit and save all models: ")
        + args.method,
        "--config",
        str(config_path),
    ]
    for path in files:
        command += ["--input", str(path)]
    result = subprocess.run(
        command + ["--", "python", "scripts/external_tabicl_readout.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    with log_path.open("x") as stream:
        stream.write(result.stdout + result.stderr)
    print(result.stdout + result.stderr, flush=True)
    result.check_returncode()
    runs = re.findall(r"run_id=(\S+)", result.stdout)
    assert len(runs) == 1
    print(json.dumps(dict(run_id=runs[0], config=str(config_path), log=str(log_path))), flush=True)


if __name__ == "__main__":
    main()
