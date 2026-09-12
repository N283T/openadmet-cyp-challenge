"""Register exact fixed-five OOF components and immutable source audit references."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--tuned", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
reg = root / "artifacts/runs/regression"
shared = root / "artifacts/runs/shared"
files = []


def add(path):
    assert path.is_file(), path
    if path not in files:
        files.append(path)
    return f"{files.index(path):03d}-{path.name}"


blend = reg / (
    "20260912T072935321335Z-1e8d2643" if args.tuned else "20260912T072501360501Z-f8c3a807"
)
audit = reg / (
    "20260912T072935933157Z-6e1a3bce" if args.tuned else "20260912T072502729280Z-d8928353"
)
cfg = dict(
    final_oof=add(blend / "outputs/external-oof.parquet"),
    blend_config=add(
        root
        / "configs/external-ensemble"
        / (
            "caruana_medium_scores_tuned_ft_cheme_v1.json"
            if args.tuned
            else "caruana_medium_scores_ft_cheme_equal_v1.json"
        )
    ),
    blend_manifest=add(blend / "manifest.json"),
    prior_macro=add(audit / "outputs/macro.csv"),
    measurements=add(shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"),
    assignments=add(reg / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"),
    scorer=add(root / "artifacts/metric-reference/custom_scoring_functions.py"),
    components=[],
    references=[],
    ft_membership=[],
)
for run, name in [
    ("20260912T054406239560Z-96fce32b", "Original nested Caruana"),
    ("20260912T060059705105Z-2f0fbfbd", "Jeremy medium frozen embedding TabPFN"),
    ("20260912T072303592078Z-07f2749a", "Jeremy upstream FFN scores plus nested log2fc TabPFN"),
    (
        "20260912T061629455493Z-730c8bea",
        "Jeremy medium fine-tuned embedding plus nested log2fc LGBM",
    ),
    ("20260912T071619253798Z-73a6fb7e", "Jeremy CheMeleon PCA256 plus nested log2fc TabICL"),
]:
    if args.tuned and run.endswith("730c8bea"):
        matches = list(reg.glob("*-e6b6cbcb"))
        assert len(matches) == 1
        run = matches[0].name
        name += " (inner-search tuned)"
    d = reg / run
    item = dict(name=name, run_id=run, oof=add(d / "outputs/external-oof.parquet"))
    for f in ["manifest", "intent", "result", "config"]:
        item[f] = add(d / (f + ".json"))
    cfg["components"].append(item)
for run in [
    "20260912T060127278065Z-39e5a986",
    "20260912T060426566774Z-7604ffeb",
    "20260912T055819206227Z-c060d197",
    "20260912T061231617816Z-299439e8",
    "20260912T052549384296Z-bac5fb3a",
]:
    d = shared / run
    cfg["references"].append(
        dict(
            run_id=run,
            manifest=add(d / "manifest.json"),
            result=add(d / "result.json"),
            summary=add(d / "outputs/summary.json"),
        )
    )
d = reg / "20260912T060538153709Z-8d434eee"
for fold in range(5):
    cfg["ft_membership"].append(add(d / f"outputs/fold{fold}/membership.parquet"))
cfg["references"].append(
    dict(
        run_id=d.name,
        manifest=add(d / "manifest.json"),
        result=add(d / "result.json"),
        summary=add(d / "outputs/summary.json"),
    )
)
cfg["nested_membership_audit"] = add(
    shared / "20260912T052549384296Z-bac5fb3a/outputs/membership-audit.json"
)
cfg["public_inventory_report"] = add(
    root / "docs/analysis/2026-09-12-jeremy-pretraining-public-inventory.md"
)
folder = root / "configs/final-regression"
folder.mkdir(exist_ok=True)
config = folder / (
    "independent-oof-audit-tuned-v1.json" if args.tuned else "independent-oof-audit-v1.json"
)
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
    "Independent final fixed-five OOF reconstruction, official scoring and provenance package",
    "--config",
    str(config),
]
for path in files:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/finalize_regression_audit.py"], cwd=root, check=True
)
