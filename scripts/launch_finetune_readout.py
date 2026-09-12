"""Evaluate an existing completed fold-specific FT run without retraining its encoder."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit-run", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--backend", choices=["lgbm", "tabpfn", "catboost"], required=True)
    parser.add_argument("--nested", type=Path)
    parser.add_argument("--pca", type=int, default=0)
    args = parser.parse_args()
    assert re.fullmatch(r"[a-z0-9_]+", args.tag)
    assert re.fullmatch(r"\d{8}T\d{12}Z-[a-f0-9]{8}", args.fit_run)
    assert args.pca >= 0
    root = Path(__file__).resolve().parents[1]
    fit_dir = root / f"artifacts/runs/regression/{args.fit_run}"
    assert (fit_dir / "result.json").is_file(), "FT run has not finished"
    result = json.loads((fit_dir / "result.json").read_text())
    assert result.get("status") == "completed" and result.get("exit_code") == 0, result
    fit_output = fit_dir / "outputs"
    fit_summary = json.loads((fit_output / "summary.json").read_text())
    assert not fit_summary.get("preparation_only", False)
    assert fit_summary["oof_rows"] == 6384
    split_run = fit_summary["config"]["split_run"]
    assert split_run == "20260911T104703655337Z-a3923255"
    data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    assignments = root / f"artifacts/runs/regression/{split_run}/outputs/parent-assignments.parquet"
    inputs = [fit_output / f"fold{fold}/external-embeddings.parquet" for fold in range(5)]
    config = dict(
        method=args.tag,
        backend=args.backend,
        feature_prefix="emb_",
        pca=args.pca,
        split_run=split_run,
        expected_oof_rows=6384,
        finetune_run=args.fit_run,
        feature_files={str(fold): f"{fold:03d}-external-embeddings.parquet" for fold in range(5)},
    )
    inputs.extend([data / "measurements.parquet", assignments])
    if args.nested:
        inputs.append(args.nested.resolve())
        config["nested_file"] = args.nested.name
    if args.backend == "tabpfn":
        weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
        config["weights"] = f"{len(inputs):03d}-{weights.name}"
        inputs.append(weights)
    baseline = (
        root
        / "artifacts/runs/regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet"
    )
    scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"
    assert all(path.is_file() for path in inputs + [baseline, scorer])
    configs, logs = root / "configs/external-finetune", root / "artifacts/launches"
    configs.mkdir(exist_ok=True)
    logs.mkdir(exist_ok=True)
    manifest_path = logs / (args.tag + "_readout_manifest.json")
    assert not manifest_path.exists()
    for name in [args.tag, args.tag + "_audit"]:
        assert not (configs / (name + ".json")).exists()
        assert not (logs / (name + ".log")).exists()

    def launch(name, configuration, declared_inputs, script, purpose):
        config_path = configs / (name + ".json")
        with config_path.open("x") as stream:
            json.dump(configuration, stream, indent=2)
        command = [
            "pixi",
            "run",
            "cyp",
            "run",
            "--track",
            "regression",
            "--purpose",
            purpose,
            "--config",
            str(config_path),
        ]
        for path in declared_inputs:
            command += ["--input", str(path)]
        process = subprocess.run(
            command + ["--", "python", script], cwd=root, capture_output=True, text=True
        )
        with (logs / (name + ".log")).open("x") as stream:
            stream.write(process.stdout + process.stderr)
        print(process.stdout, flush=True)
        process.check_returncode()
        matches = re.findall(r"run_id=(\S+)", process.stdout)
        assert len(matches) == 1
        return matches[0]

    run = launch(
        args.tag,
        config,
        inputs,
        "scripts/readout_external_finetune.py",
        "Readout of completed fold-specific FT embeddings: " + args.tag,
    )
    audit = launch(
        args.tag + "_audit",
        dict(
            model_run=run,
            database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
        ),
        [root / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet", baseline, scorer],
        "scripts/audit_external_readout.py",
        "Matched OOF and equal blend audit: " + args.tag,
    )
    manifest = dict(
        finetune_run=args.fit_run,
        readout_run=run,
        audit_run=audit,
        backend=args.backend,
        method=args.tag,
    )
    with manifest_path.open("x") as stream:
        json.dump(manifest, stream, indent=2)
    print(json.dumps(manifest), flush=True)


if __name__ == "__main__":
    main()
