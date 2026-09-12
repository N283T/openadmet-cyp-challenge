"""Declare artifacts for fixed seed and scope means, without score selection."""

import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    inputs = []

    def add(path):
        assert path.is_file(), path
        name = f"{len(inputs):03d}-{path.name}"
        inputs.append(path)
        return name

    split = "20260911T104703655337Z-a3923255"
    global_run = "20260912T094021204217Z-84c28266"
    cfg = dict(
        split_run=split,
        global_run=global_run,
        predictions={},
        means={},
        memberships={},
        histories={},
    )
    cfg["global_file"] = add(
        root / f"artifacts/runs/regression/{global_run}/outputs/global-oof.parquet"
    )
    cfg["assignments"] = add(
        root / f"artifacts/runs/regression/{split}/outputs/parent-assignments.parquet"
    )
    tags = {
        "seed42": "d6_aux",
        "seed17": "d6_aux_seed17",
        "seed89": "d6_aux_seed89",
        "full": "full",
        "direct": "d6_direct",
        "freshscratch": "d6_direct_fresh_scratchrates",
    }
    for name, tag in tags.items():
        manifest = json.loads(
            (root / f"artifacts/launches/d6_specialist_{tag}_v1_manifest.json").read_text()
        )
        fit = root / f"artifacts/runs/regression/{manifest['fit_run']}/outputs"
        ro = root / f"artifacts/runs/regression/{manifest['readout_run']}/outputs"
        cfg["predictions"][name + "_direct"] = add(fit / "finetune-oof.parquet")
        cfg["predictions"][name + "_lgbm"] = add(ro / "external-oof.parquet")
        if name.startswith("seed"):
            seed = name.removeprefix("seed")
            cfg["memberships"][seed] = [
                add(fit / f"fold{fold}/membership.parquet") for fold in range(5)
            ]
            cfg["histories"][seed] = [add(fit / f"fold{fold}/history.jsonl") for fold in range(5)]
    for backend in ["direct", "lgbm"]:
        cfg["means"]["d6_seed_mean_" + backend] = [
            "seed" + str(seed) + "_" + backend for seed in [42, 17, 89]
        ]
        cfg["means"]["d6_scope_mean_" + backend] = [
            name + "_" + backend for name in ["full", "seed42", "direct", "freshscratch"]
        ]
    add(root / "artifacts/metric-reference/custom_scoring_functions.py")
    config = root / "configs/d6-specialist/seed-summary-v1.json"
    with config.open("x") as stream:
        json.dump(cfg, stream, indent=2)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Prespecified3seed and4scope specialist means, membership audits and fixedglobalblends",
        "--config",
        str(config),
    ]
    for item in inputs:
        command += ["--input", str(item)]
    result = subprocess.run(
        command + ["--", "python", "scripts/d6_specialist_seed_summary.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    with (root / "artifacts/launches/d6_specialist_seed_summary_v1.log").open("x") as stream:
        stream.write(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()


if __name__ == "__main__":
    main()
