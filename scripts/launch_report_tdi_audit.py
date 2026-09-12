"""Launch the immutable paired TDI experiment comparison after both arms finish."""

import json
import re
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    manifests = {
        name: json.loads(
            (root / f"artifacts/launches/report_tdi_{name}_v1_ft_manifest.json").read_text()
        )
        for name in ["arm", "arm_delta"]
    }
    baseline_run = "20260912T083627802028Z-ae562a59"
    inputs = [
        root / f"artifacts/runs/regression/{baseline_run}/outputs/focus_2d6_combined_50pct.parquet"
    ]
    cfg = dict(
        baseline_run=baseline_run,
        baseline_file="000-" + inputs[0].name,
        prediction_files={},
        membership_files={},
    )

    def add(path):
        name = f"{len(inputs):03d}-{path.name}"
        inputs.append(path)
        return name

    control_run = "20260912T072355846776Z-2e678154"
    control_readout = "20260912T092824781073Z-fef7788e"
    models = {"control": (control_run, control_readout)}
    models.update(
        {
            name: (m["finetune_run"], m["readouts"]["lgbm"]["model_run"])
            for name, m in manifests.items()
        }
    )
    for name, (fit, readout) in models.items():
        fit_out = root / f"artifacts/runs/regression/{fit}/outputs"
        cfg["prediction_files"][name + "_direct"] = add(fit_out / "finetune-oof.parquet")
        cfg["prediction_files"][name + "_lgbm"] = add(
            root / f"artifacts/runs/regression/{readout}/outputs/external-oof.parquet"
        )
        cfg["membership_files"][name] = [
            add(fit_out / f"fold{fold}/membership.parquet") for fold in range(5)
        ]
    add(root / "artifacts/metric-reference/custom_scoring_functions.py")
    config = root / "configs/external-finetune/report_tdi_comparison_v2.json"
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
        "Compare continuous TDI auxiliary arms and fixed25/50percent ensemble additions; "
        "audit all outer memberships",
        "--config",
        str(config),
    ]
    for path in inputs:
        command += ["--input", str(path)]
    proc = subprocess.run(
        command + ["--", "python", "scripts/report_tdi_audit.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    log = root / "artifacts/launches/report_tdi_comparison_v2.log"
    with log.open("x") as stream:
        stream.write(proc.stdout + proc.stderr)
    print(proc.stdout, flush=True)
    proc.check_returncode()
    run = re.findall(r"run_id=(\S+)", proc.stdout)
    assert len(run) == 1
    with (root / "artifacts/launches/report_tdi_comparison_v2_manifest.json").open("x") as stream:
        json.dump(dict(audit_run=run[0], models=models), stream, indent=2)


if __name__ == "__main__":
    main()
