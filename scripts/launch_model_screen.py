"""Persist configurations and launch a bounded feature/direct-FT comparison."""

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/shared"
STACK = [
    "20260910T121031603901Z-fbc2b2ea",
    "20260910T121149831493Z-a62f51ef",
    "20260910T121402665787Z-08f7afb5",
    "20260910T121548828209Z-abe37170",
    "20260910T121624900089Z-271684c2",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("variant", choices=["mordred", "stack_mordred", "direct_ft"])
    args = parser.parse_args()
    for fold in range(5):
        cfg = json.loads((ROOT / f"configs/readout/tabpfn-embedding-fold{fold}.json").read_text())
        cfg.pop("encoder_run")
        cfg.update(variant=args.variant, feature_run="20260909T132918829974Z-bc81af58")
        files = [
            RUNS / cfg["dataset_run"] / "outputs/measurements.parquet",
            RUNS / cfg["split_run"] / "outputs/assignments.parquet",
        ]
        if args.variant == "direct_ft":
            cfg.update(
                max_epochs=100,
                patience=15,
                inner_validation_fraction=0.1,
                batch_size=64,
                learning_rate=1e-4,
                ffn_hidden_dim=256,
                ffn_dropout=0.1,
                initialization="chemeleon",
                selection="inner_loss",
                classification_threshold=0.5,
            )
            files += [
                RUNS / cfg["dataset_run"] / "outputs/compounds.parquet",
                RUNS / "20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt",
            ]
            script = "scripts/direct_ft.py"
        else:
            files += sorted((RUNS / cfg["feature_run"] / "outputs").glob("mordred-*.parquet"))
            files += [
                Path.home() / f".cache/tabpfn/tabpfn-v2.6-{k}-v2.6_default.ckpt"
                for k in ["regressor", "classifier"]
            ]
            if args.variant == "stack_mordred":
                cfg["stack_run"] = STACK[fold]
                files += sorted((RUNS / STACK[fold] / "outputs").glob("*-features.parquet"))
                files += [RUNS / STACK[fold] / "outputs/summary.json"]
            script = "scripts/descriptor_readout.py"
        path = ROOT / f"configs/screen/{args.variant}-fold{fold}.json"
        path.parent.mkdir(exist_ok=True)
        assert not path.exists()
        path.write_text(json.dumps(cfg, indent=2) + "\n")
        cmd = [
            "pixi",
            "run",
            "cyp",
            "run",
            "--track",
            "shared",
            "--purpose",
            f"{args.variant} UMAP outer fold {fold}; fixed evaluation protocol",
            "--config",
            str(path),
        ]
        for p in files:
            cmd += ["--input", str(p)]
        subprocess.run(cmd + ["--", "python", script], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
