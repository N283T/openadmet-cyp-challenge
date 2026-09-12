"""Record official standalone audit and fixed three-way FT HPO comparison."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = "20260912T072056445264Z-e6b6cbcb"
outputs = root / "artifacts/runs/regression" / run / "outputs"
scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"


def launch(tag, cfg, inputs, script):
    config = root / "configs/external-ft-tuning" / (tag + ".json")
    with config.open("x") as f:
        json.dump(cfg, f, indent=2)
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        tag.replace("_", " "),
        "--config",
        str(config),
    ]
    for p in inputs:
        cmd += ["--input", str(p)]
    result = subprocess.run(
        cmd + ["--", "python", script], cwd=root, capture_output=True, text=True
    )
    print(result.stdout, flush=True)
    print(result.stderr, flush=True)
    result.check_returncode()
    return re.search(r"run_id=(\S+)", result.stdout).group(1)


launch(
    "jeremy_ft_logfc_rustuna_v1_audit",
    dict(
        model_run=run,
        database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
    ),
    [
        outputs / "external-oof.parquet",
        root
        / "artifacts/runs/regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet",
        scorer,
    ],
    "scripts/audit_external_readout.py",
)
inputs = []


def add(path):
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


cfg = dict(
    fit_run=run,
    models={},
    scorer=add(scorer),
    membership=add(outputs / "membership.parquet"),
    summary=add(outputs / "summary.json"),
    encoder_exclusion_audit=add(outputs / "encoder-exclusion-audit.json"),
    trial_files=[],
)
for name, r in {
    "caruana": "20260912T054406239560Z-96fce32b",
    "frozen_pfn": "20260912T060059705105Z-2f0fbfbd",
    "fixed_lgbm": "20260912T061629455493Z-730c8bea",
    "tuned_lgbm": run,
}.items():
    cfg["models"][name] = add(root / f"artifacts/runs/regression/{r}/outputs/external-oof.parquet")
cfg["reference_membership"] = add(
    root / "artifacts/runs/regression/20260912T054406239560Z-96fce32b/outputs/membership.parquet"
)
for p in sorted(outputs.glob("fold-*/*/trials.jsonl")):
    cfg["trial_files"].append(
        dict(file=add(p), fold=int(p.parent.parent.name.split("-")[1]), endpoint=p.parent.name)
    )
assert len(cfg["trial_files"]) == 20
r = launch(
    "jeremy_ft_logfc_rustuna_v1_report", cfg, inputs, "scripts/audit_tuned_finetuned_readout.py"
)
report = root / f"artifacts/runs/regression/{r}/outputs/report.md"
with (root / "docs/analysis/2026-09-12-jeremy-ft-logfc-rustuna.md").open("x") as f:
    f.write(report.read_text())
