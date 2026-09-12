"""Launch fixed-scope specialist comparison after all five fits finish."""

import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    run = "20260912T094021204217Z-84c28266"
    inputs = [root / f"artifacts/runs/regression/{run}/outputs/global-oof.parquet"]
    cfg = dict(
        global_run=run,
        global_file="000-global-oof.parquet",
        models={},
        split_run="20260911T104703655337Z-a3923255",
    )

    def add(path):
        assert path.is_file()
        name = f"{len(inputs):03d}-{path.name}"
        inputs.append(path)
        return name

    for suffix in [
        "full",
        "d6_aux",
        "d6_direct",
        "d6_direct_fresh",
        "d6_direct_fresh_scratchrates",
    ]:
        tag = "d6_specialist_" + suffix + "_v1"
        manifest = json.loads((root / "artifacts/launches" / (tag + "_manifest.json")).read_text())
        fit = root / f"artifacts/runs/regression/{manifest['fit_run']}/outputs"
        readout = root / f"artifacts/runs/regression/{manifest['readout_run']}/outputs"
        cfg["models"][suffix] = dict(
            scope=manifest["scope"],
            fit_run=manifest["fit_run"],
            readout_run=manifest["readout_run"],
            summary=add(fit / "summary.json"),
            direct=add(fit / "finetune-oof.parquet"),
            lgbm=add(readout / "external-oof.parquet"),
            memberships=[add(fit / f"fold{fold}/membership.parquet") for fold in range(5)],
            histories={str(fold): add(fit / f"fold{fold}/history.jsonl") for fold in range(5)},
        )
    add(root / "artifacts/metric-reference/custom_scoring_functions.py")
    add(root / f"artifacts/runs/regression/{cfg['split_run']}/outputs/parent-assignments.parquet")
    config = root / "configs/d6-specialist/comparison-v2.json"
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
        "Audit specialist scopes and fixed additions to latest global OOF",
        "--config",
        str(config),
    ]
    for item in inputs:
        command += ["--input", str(item)]
    result = subprocess.run(
        command + ["--", "python", "scripts/d6_specialist_audit.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    with (root / "artifacts/launches/d6_specialist_comparison_v2.log").open("x") as stream:
        stream.write(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()


if __name__ == "__main__":
    main()
