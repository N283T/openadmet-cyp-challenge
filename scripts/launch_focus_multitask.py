"""Launch immutable external600 multitask fine-tuning, readouts and matched audits."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--logfc-weight", type=float, default=0.25)
    parser.add_argument(
        "--backends", nargs="+", choices=["lgbm", "tabpfn"], default=["lgbm", "tabpfn"]
    )
    parser.add_argument("--epochs", type=int, default=60)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--pca", type=int, default=0)
    parser.add_argument("--extra-targets", type=Path)
    parser.add_argument("--extra-weight", type=float, default=1.0)
    parser.add_argument("--nested", type=Path)
    parser.add_argument("--training-loss", choices=["mse", "mae", "pic50_interval"], default="mse")
    parser.add_argument("--auxiliary-loss", choices=["mse", "mae"], default="mse")
    parser.add_argument(
        "--checkpoint-selection",
        choices=["multitask_mse", "pic50_st_rae"],
        default="multitask_mse",
    )
    args = parser.parse_args()
    assert re.fullmatch(r"[a-z0-9_]+", args.tag)
    assert len(args.backends) == len(set(args.backends))
    root = Path(__file__).resolve().parents[1]
    checkpoint = args.checkpoint.resolve()
    data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    split_run = "20260911T104703655337Z-a3923255"
    assignments = root / f"artifacts/runs/regression/{split_run}/outputs/parent-assignments.parquet"
    targets = root / "artifacts/runs/shared/20260910T100337258849Z-927eac1c/outputs/targets.parquet"
    baseline = (
        root
        / "artifacts/runs/regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet"
    )
    scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"
    weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
    required = [
        checkpoint,
        data / "compounds.parquet",
        data / "measurements.parquet",
        assignments,
        targets,
        baseline,
        scorer,
    ]
    if "tabpfn" in args.backends:
        required.append(weights)
    if args.extra_targets:
        required.append(args.extra_targets.resolve())
    if args.nested:
        required.append(args.nested.resolve())
    assert all(path.is_file() for path in required), required
    config_dir, log_dir = root / "configs/external-finetune", root / "artifacts/launches"
    config_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)
    names = [args.tag + "_ft", args.tag + "_ft_direct_audit"]
    for backend in args.backends:
        names.extend([args.tag + "_ft_" + backend, args.tag + "_ft_" + backend + "_audit"])
    manifest_path = log_dir / (args.tag + "_ft_manifest.json")
    assert not manifest_path.exists()
    for name in names:
        assert not (config_dir / (name + ".json")).exists()
        assert not (log_dir / (name + ".log")).exists()

    def launch(name, config, inputs, script, purpose):
        assert all(path.is_file() for path in inputs)
        config_path = config_dir / (name + ".json")
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
            purpose,
            "--config",
            str(config_path),
        ]
        for path in inputs:
            command += ["--input", str(path)]
        result = subprocess.run(
            command + ["--", "python", script], cwd=root, capture_output=True, text=True
        )
        with (log_dir / (name + ".log")).open("x") as stream:
            stream.write(result.stdout + result.stderr)
        print(result.stdout, flush=True)
        result.check_returncode()
        runs = re.findall(r"run_id=(\S+)", result.stdout)
        assert len(runs) == 1
        return runs[0]

    def audit(name, run, filename):
        return launch(
            name + "_audit",
            dict(
                model_run=run,
                prediction_file=filename,
                database=dict(
                    host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"
                ),
            ),
            [root / f"artifacts/runs/regression/{run}/outputs/{filename}", baseline, scorer],
            "scripts/audit_external_readout.py",
            "Matched official scores and prespecified equal blends: " + name,
        )

    fit_config = dict(
        training_loss=args.training_loss,
        auxiliary_loss=args.auxiliary_loss,
        checkpoint_selection=args.checkpoint_selection,
        seed=42,
        encoder_lr=1e-5,
        head_lr=1e-4,
        weight_decay=1e-4,
        max_epochs=args.epochs,
        patience=args.patience,
        batch_size=64,
        inference_batch_size=128,
        dropout=0.1,
        inner_fraction=0.1,
        inner_group_column="group",
        device="cuda",
        num_threads=4,
        pIC50_weight=2.0,
        pic50_endpoint_weights={"CYP1A2": 2.0, "CYP2C9": 0.5, "CYP2D6": 2.0, "CYP3A4": 0.5},
        selection_targets=["CYP1A2", "CYP2D6"],
        log2fc_weight=args.logfc_weight,
        checkpoint_file=checkpoint.name,
        split_run=split_run,
        expected_oof_rows=6384,
        method=args.tag + "_ft_direct",
        source_checkpoint=str(checkpoint),
        outer_policy="Exclude every label of outer nonanchor compounds; anchors train",
    )
    fit_inputs = [
        checkpoint,
        data / "compounds.parquet",
        data / "measurements.parquet",
        assignments,
        targets,
    ]
    if args.extra_targets:
        fit_config["auxiliary_file"] = "004-targets.parquet"
        extra_filename = f"{len(fit_inputs):03d}-{args.extra_targets.name}"
        fit_inputs.append(args.extra_targets.resolve())
        fit_config.update(extra_targets_file=extra_filename, extra_weight=args.extra_weight)
    if args.checkpoint_selection == "pic50_st_rae" or args.training_loss == "pic50_interval":
        fit_inputs.append(scorer)
    fit_run = launch(
        args.tag + "_ft",
        fit_config,
        fit_inputs,
        "scripts/focus_multitask.py",
        "External600 -> pIC50/log2fc and configured auxiliary heads; outer compound exclusion: "
        + args.tag,
    )
    manifest = dict(
        tag=args.tag,
        checkpoint=str(checkpoint),
        finetune_run=fit_run,
        direct_audit=audit(args.tag + "_ft_direct", fit_run, "finetune-oof.parquet"),
        readouts={},
    )
    fit_output = root / f"artifacts/runs/regression/{fit_run}/outputs"
    for backend in args.backends:
        name = args.tag + "_ft_" + backend
        inputs = [fit_output / f"fold{fold}/external-embeddings.parquet" for fold in range(5)]
        config = dict(
            method=name,
            backend=backend,
            feature_prefix="emb_",
            pca=args.pca,
            split_run=split_run,
            expected_oof_rows=6384,
            finetune_run=fit_run,
            feature_files={
                str(fold): f"{fold:03d}-external-embeddings.parquet" for fold in range(5)
            },
        )
        inputs.extend([data / "measurements.parquet", assignments])
        if args.nested:
            inputs.append(args.nested.resolve())
            config["nested_file"] = args.nested.name
        if backend == "tabpfn":
            config["weights"] = f"{len(inputs):03d}-{weights.name}"
            inputs.append(weights)
        readout_run = launch(
            name,
            config,
            inputs,
            "scripts/readout_external_finetune.py",
            "Fold-specific adapted external600 readout: " + name,
        )
        manifest["readouts"][backend] = dict(
            model_run=readout_run, audit_run=audit(name, readout_run, "external-oof.parquet")
        )
    with manifest_path.open("x") as stream:
        json.dump(manifest, stream, indent=2)
    print(json.dumps(manifest), flush=True)


if __name__ == "__main__":
    main()
