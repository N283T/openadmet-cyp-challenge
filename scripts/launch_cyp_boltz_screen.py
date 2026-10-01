"""Launch reproducible full-five-fold trunk pooling/readout comparisons."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("target", choices=["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"])
    parser.add_argument("pool_run")
    parser.add_argument("--track", choices=["regression", "classification"], default="regression")
    parser.add_argument("--protocol", choices=["umap", "parent-ensemble"], default="umap")
    args = parser.parse_args()
    if args.track == "classification" and args.target not in ["cyp2d6", "cyp3a4"]:
        raise ValueError("No TDI endpoint for this CYP")
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/shared"
    weights = runs / "20260910T122343854890Z-3a03f908/inputs"
    kind = "classifier" if args.track == "classification" else "regressor"
    files = {
        "trunk.parquet": runs / args.pool_run / "outputs/features.parquet",
        "pool-summary.json": runs / args.pool_run / "outputs/summary.json",
        "pool-config.json": runs / args.pool_run / "config.json",
        "assignments.parquet": runs / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
        "measurements.parquet": runs
        / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        "tabpfn.ckpt": next(weights.glob(f"*tabpfn-v2.6-{kind}-v2.6_default.ckpt")),
    }
    for i, path in enumerate(
        sorted((runs / "20260910T012630927160Z-e9576045/outputs").glob("chemeleon-*.parquet"))
    ):
        files[f"cheme-{i:04d}.parquet"] = path
    for i, path in enumerate(
        sorted((runs / "20260909T132918829974Z-bc81af58/outputs").glob("mordred-*.parquet"))
    ):
        files[f"mordred-{i:04d}.parquet"] = path
    suffix = "v1" if args.protocol == "umap" else "parent-ensemble-v2"
    config = root / "configs/boltz-trunk" / f"screen-{args.target}-{args.track}-{suffix}.json"
    if not config.exists():
        endpoint = args.target.upper() + (
            "_is_TDI" if args.track == "classification" else "_pIC50_direct_inhibition"
        )
        config.write_text(
            json.dumps(
                {
                    "target": args.target,
                    "track": args.track,
                    "endpoint": endpoint,
                    "split_id": "20260909T135130268493Z-edba6fc0:umap",
                    "folds": [0, 1, 2, 3, 4],
                    "seed": 20261001,
                    "threads": 4,
                    "tabpfn_version": "2.6",
                    "tabpfn_estimators": 8,
                    "feature_selection": "LGBM training-fold gain top500, no outer-label use",
                    "views": [
                        "ligand_only",
                        "pxr_protein",
                        "pxr_core",
                        "pxr_extended",
                        "core_mixed_hem",
                        "core_separate_hem",
                        "extended_mixed_hem",
                        "extended_separate_hem",
                        "core_separate_fe",
                        "core_bidirectional_hem",
                        "core_separate_hem_s",
                        "core_distribution_hem",
                    ],
                    "feature_modes": ["trunk", "trunk_molecular"],
                    "molecular_control": "frozen CheMeleon2048 + Mordred1613",
                    "readouts": ["lgbm_binary", "tabpfn_top500"]
                    if args.track == "classification"
                    else ["lgbm_l1", "lgbm_l2", "tabpfn_top500"],
                    "blend_weights": [0.1, 0.25, 0.5],
                    "selection": "Exploratory fixed grid; no candidate promotion",
                    "classification_scope": (
                        "Incomplete all-molecule pool: diagnostic only; no call-rate selection"
                    ),
                },
                indent=2,
            )
            + "\n"
        )
    if args.protocol == "parent-ensemble":
        if args.track != "regression":
            raise ValueError("TDI uses existing matched UMAP runs, without refitting")
        regression = root / "artifacts/runs/regression"
        controls = {
            "cyp1a2": (
                "20260918T034159381060Z-46314b50/outputs/1A2/add_dense_and_fast_pool/oof.parquet",
                {"method": "equal_all"},
                1372,
            ),
            "cyp2c9": (
                "20260917T123601175420Z-e1560701/outputs/oof.parquet",
                {"method": "equal_subset"},
                1255,
            ),
            "cyp2d6": (
                "20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
                {
                    "split_id": "parent",
                    "robust_delta": 0.5,
                    "constraint_weight": 0.35,
                    "placement": "fold_moment_match",
                    "blend_weight": 0.5,
                },
                1457,
            ),
            "cyp3a4": (
                "20260917T130718808881Z-dfa02730/outputs/oof.parquet",
                {"method": "caruana"},
                2300,
            ),
        }
        source, selection, count = controls[args.target]
        files["control-oof.parquet"] = regression / source
        files["assignments.parquet"] = (
            regression / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
        )
        cfg = json.loads(config.read_text())
        cfg.update(
            {
                "cohort_from_control": True,
                "control_filter": selection,
                "expected_labelled_rows": count,
                "assignment_strategy": "series",
                "split_id": "20260911T104703655337Z-a3923255:series",
                "selection": (
                    "Whole-member placement in submitted ensemble; fixed parent cohort and folds"
                ),
                "control_source": str(files["control-oof.parquet"].relative_to(root)),
            }
        )
        config.write_text(json.dumps(cfg, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        args.track,
        "--purpose",
        (
            f"Compare 12 PXR-style CYP/HEM pooling views with LGBM and TabPFN "
            f"for {args.target} {args.track}"
        ),
        "--config",
        str(config),
    ]
    for alias, path in files.items():
        command += ["--input", f"{alias}={path}"]
    subprocess.run(
        command + ["--", "python", "scripts/screen_cyp_boltz_trunk.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
