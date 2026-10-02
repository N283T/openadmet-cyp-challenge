"""Reuse completed trunk TDI fits for exact submitted whole-member placement."""

import json
import sys
from pathlib import Path

from launch_cyp_boltz_incumbent_suite import recorded


def main():
    root = Path(__file__).resolve().parents[1]
    directory = Path(sys.argv[1])
    runs = root / "artifacts/runs/classification"
    d6 = runs / "20260925T130752164378Z-d40a0d9c"
    a4 = runs / "20260915T125001009276Z-d5adff22"
    inputs = {
        "depth-oof.parquet": root
        / "artifacts/runs/shared/20260925T065626160287Z-8017d516/outputs/oof.parquet",
        "d6-current.parquet": d6 / "outputs/crossfit-mean-block-oof.parquet",
        "d6-lgbm.parquet": d6 / "inputs/lgbm-oof.parquet",
        "d6-pfn.parquet": d6 / "inputs/pfn-oof.parquet",
        "d6-graph.parquet": d6 / "inputs/graph-oof.parquet",
        "d6-fixed.parquet": d6 / "inputs/fixed-oof.parquet",
        "d6-raw600.parquet": d6 / "inputs/raw600-oof.parquet",
        "d6-selections.parquet": d6 / "inputs/selections.parquet",
        "a4-oof.parquet": a4 / "outputs/oof.parquet",
        "a4-member.parquet": a4 / "inputs/000-member-oof.parquet",
        "a4-raw.parquet": a4 / "inputs/004-outer-oof.parquet",
        "a4-aux.parquet": a4 / "inputs/002-oof.parquet",
        "a4-selections.parquet": a4 / "outputs/selections.parquet",
        "d6-trunk.parquet": runs / "20261001T023655009132Z-91096754/outputs/predictions.parquet",
        "a4-trunk.parquet": runs / "20261001T043809889682Z-cd81b014/outputs/predictions.parquet",
    }
    cfg = json.loads(
        (root / "configs/screen/classification-depth-whole-members-v2.json").read_text()
    )
    cfg.update(
        {
            "depth_sets": [],
            "backends": {},
            "trunk_predictions": {"d6": "d6-trunk.parquet", "a4": "a4-trunk.parquet"},
            "scope": (
                "Exact submitted UMAP OOF whole-member addition/replacement; "
                "AP-selected placement excludes held fold; labelled diagnostic only; "
                "no deployment threshold selection"
            ),
        }
    )
    identifier = recorded(
        directory,
        "tdi-incumbent-placement",
        "classification",
        cfg,
        inputs,
        "evaluate_classification_depth_whole_members.py",
    )
    (directory / "tdi-placement-run.txt").write_text(identifier + "\n")
    print(f"tdi_placement_run={identifier}", flush=True)


if __name__ == "__main__":
    main()
