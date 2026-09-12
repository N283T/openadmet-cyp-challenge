"""Recorded scope-matched CYP2D6 neural specialist and readout experiments."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=["full", "d6_aux", "d6_direct"], required=True)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--scratch-rates", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    assert not args.scratch_rates or args.fresh
    root = Path(__file__).resolve().parents[1]
    configs = root / "configs/d6-specialist"
    configs.mkdir(exist_ok=True)
    data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    split = "20260911T104703655337Z-a3923255"
    assignments = root / f"artifacts/runs/regression/{split}/outputs/parent-assignments.parquet"
    inputs = [
        root / "artifacts/runs/shared/20260912T060426566774Z-7604ffeb/outputs/best.pt",
        data / "compounds.parquet",
        data / "measurements.parquet",
        assignments,
        root
        / "artifacts/runs/shared/20260912T093252891365Z-271e3020/outputs"
        / "observation-targets.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ]
    cfg = dict(
        scope=args.scope,
        fresh_initialization=args.fresh,
        seed=args.seed,
        encoder_lr=1e-3 if args.scratch_rates else 1e-5,
        head_lr=1e-3 if args.scratch_rates else 1e-4,
        weight_decay=1e-4,
        max_epochs=160,
        patience=20,
        batch_size=64,
        inference_batch_size=128,
        dropout=0.1,
        inner_fraction=0.1,
        inner_group_column="group",
        device="cuda",
        num_threads=2,
        pIC50_weight=2.0,
        log2fc_weight=1.0,
        aux_kind_weights=dict(tdi=0.5, emax=0.25, logfc=1.0),
        selection_targets=["CYP2D6"],
        checkpoint_file="best.pt",
        training_loss="pic50_interval",
        auxiliary_loss="mse",
        checkpoint_selection="pic50_st_rae",
        split_run=split,
        expected_oof_rows=6384 if args.scope == "full" else 1457,
        method="d6_specialist_" + args.scope + "_direct",
    )
    if args.seed != 42:
        cfg["inner_split_seed"] = 42

    def launch(tag, config, declared, script):
        path = configs / (tag + ".json")
        with path.open("x") as stream:
            json.dump(config, stream, indent=2)
        command = [
            "pixi",
            "run",
            "cyp",
            "run",
            "--track",
            "regression",
            "--purpose",
            "Matched CYP2D6 specialist scope: " + tag,
            "--config",
            str(path),
        ]
        for item in declared:
            assert item.is_file()
            command += ["--input", str(item)]
        p = subprocess.run(
            command + ["--", "python", script], cwd=root, capture_output=True, text=True
        )
        with (root / "artifacts/launches" / (tag + ".log")).open("x") as stream:
            stream.write(p.stdout + p.stderr)
        print(p.stdout, flush=True)
        p.check_returncode()
        run = re.findall(r"run_id=(\S+)", p.stdout)
        assert len(run) == 1
        return run[0]

    tag = (
        "d6_specialist_"
        + args.scope
        + ("_fresh" if args.fresh else "")
        + ("_scratchrates" if args.scratch_rates else "")
        + (f"_seed{args.seed}" if args.seed != 42 else "")
        + "_v1"
    )
    fit = launch(tag, cfg, inputs, "scripts/d6_specialist_fit.py")
    output = root / f"artifacts/runs/regression/{fit}/outputs"
    read_inputs = [output / f"fold{fold}/external-embeddings.parquet" for fold in range(5)]
    read_inputs += [
        data / "measurements.parquet",
        assignments,
        root
        / "artifacts/runs/shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    ]
    rcfg = dict(
        method=tag + "_lgbm",
        backend="lgbm",
        feature_prefix="emb_",
        pca=0,
        split_run=split,
        expected_oof_rows=1457,
        endpoint="CYP2D6_pIC50_direct_inhibition",
        finetune_run=fit,
        feature_files={str(f): f"{f:03d}-external-embeddings.parquet" for f in range(5)},
        nested_file="nested-log2fc.parquet",
    )
    readout = launch(tag + "_readout", rcfg, read_inputs, "scripts/focus_single_readout.py")
    with (root / "artifacts/launches" / (tag + "_manifest.json")).open("x") as stream:
        json.dump(dict(scope=args.scope, fit_run=fit, readout_run=readout), stream, indent=2)


if __name__ == "__main__":
    main()
