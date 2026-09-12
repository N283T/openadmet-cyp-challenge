"""Record matched official audit and independent expanded ensemble verification."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--fit-run", default="20260912T061036581812Z-288ee828")
parser.add_argument("--tag", default="jeremy_expanded_v1")
parser.add_argument("--reference-run", default="20260912T054406239560Z-96fce32b")
parser.add_argument("--classic", action="store_true")
args = parser.parse_args()
new_run = args.fit_run
old_run = "20260912T054406239560Z-96fce32b"
new = root / "artifacts/runs/regression" / new_run / "outputs"
old = root / "artifacts/runs/regression" / old_run / "outputs"
reference = root / "artifacts/runs/regression" / args.reference_run / "outputs"
baseline = (
    root / "artifacts/runs/regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet"
)
scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"


def launch(tag, cfg, paths, script):
    config = root / "configs/external-ensemble" / (tag + ".json")
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
    for p in paths:
        cmd += ["--input", str(p)]
    r = subprocess.run(cmd + ["--", "python", script], cwd=root, capture_output=True, text=True)
    print(r.stdout, flush=True)
    print(r.stderr, flush=True)
    r.check_returncode()
    return re.search(r"run_id=(\S+)", r.stdout).group(1)


launch(
    args.tag + "_audit",
    dict(
        model_run=new_run,
        database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
    ),
    [new / "external-oof.parquet", baseline, scorer],
    "scripts/audit_external_readout.py",
)
paths = []


def add(path):
    paths.append(path)
    return f"{len(paths) - 1:03d}-{path.name}"


cfg = dict(
    fit_run=new_run,
    new_oof=add(new / "external-oof.parquet"),
    old_oof=add(old / "external-oof.parquet"),
    baseline=add(baseline),
    scorer=add(scorer),
    new_membership=add(new / "membership.parquet"),
    old_membership=add(reference / "membership.parquet"),
    new_summary=add(new / "summary.json"),
    jeremy_oof=add(
        root
        / "artifacts/runs/regression/20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet"
    ),
    jeremy_ft_oof=add(
        root
        / "artifacts/runs/regression/20260912T061629455493Z-730c8bea/outputs/external-oof.parquet"
    ),
    inner_pairs=[],
    original_candidate_count=30 if args.classic else 24,
    title="Jeremy + FP/Mordred expanded nested Caruana"
    if args.classic
    else "Jeremy-expanded nested Caruana diagnostic",
)
for p in sorted(new.glob("fold-*/*/inner-predictions.parquet")):
    cfg["inner_pairs"].append(dict(new=add(p), old=add(reference / p.relative_to(new))))
assert len(cfg["inner_pairs"]) == 20
if args.classic:
    cfg["classic"] = add(
        root
        / "artifacts/runs/shared/20260912T061043373484Z-2e344d7a/outputs"
        / "fp-mordred-features.parquet"
    )
    cfg["previous_expanded_oof"] = add(reference / "external-oof.parquet")
    cfg["preprocessors"] = []
    for p in sorted(new.glob("fold-*/*/classic-*-preprocessor.npz")):
        cfg["preprocessors"].append(
            dict(
                file=add(p),
                endpoint=p.parent.name,
                outer_fold=int(p.parent.parent.name.split("-")[1]),
                role="inner" if "inner" in p.name else "outer",
            )
        )
    assert len(cfg["preprocessors"]) == 40
r = launch(args.tag + "_report", cfg, paths, "scripts/audit_jeremy_nested_ensemble.py")
report = root / f"artifacts/runs/regression/{r}/outputs/report.md"
destination = root / "docs/analysis" / ("2026-09-12-" + args.tag.replace("_", "-") + ".md")
with destination.open("x") as f:
    f.write(report.read_text())
