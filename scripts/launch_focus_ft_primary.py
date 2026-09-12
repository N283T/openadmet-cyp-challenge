"""Launch immutable fixed-parameter FT primary-feature replacement comparison."""

import itertools
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
ft = runs / "regression/20260912T060538153709Z-8d434eee/outputs"
tuned = runs / "regression/20260912T072056445264Z-e6b6cbcb/outputs"
aux = runs / "shared/20260912T082554377951Z-d272774c/outputs"
files = dict(
    assignments=runs
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    labels=runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    aux=aux / "nested-log2fc.parquet",
    old=tuned / "external-oof.parquet",
    final=runs / "regression/20260912T073049029140Z-0c26dd06/outputs/verified-oof.parquet",
    scorer=root / "artifacts/metric-reference/custom_scoring_functions.py",
)
for fold in range(5):
    files[f"emb{fold}"] = ft / f"fold{fold}/external-embeddings.parquet"
    files[f"membership{fold}"] = ft / f"fold{fold}/membership.parquet"
    for ep in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        ep += "_pIC50_direct_inhibition"
        files[f"params{fold}-{ep}"] = tuned / f"fold-{fold}/{ep}/selected.json"
for pair in itertools.combinations(range(5), 2):
    name = "".join(map(str, pair))
    files["pair" + name] = aux / f"pair{name}-predictions.parquet"
    for head in range(4):
        files[f"pair{name}-head{head}"] = aux / f"pair{name}-head{head}-membership.parquet"
assert all(p.is_file() for p in files.values())
cfg = dict(
    method="focus_ft_primary_fixedparams_v1",
    split_run="20260911T104703655337Z-a3923255",
    inputs={key: f"{i:03d}-{p.name}" for i, (key, p) in enumerate(files.items())},
    retune=False,
    fixed_additions=[0.25, 0.5],
    fixed_component_replacement=True,
)
p = root / "configs/focus-primary/ft-fixedparams-v1.json"
with p.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Replace primary features in selected FT600 LGBM component; "
        "fixed parameters, exact exclusions and 25/50 blends"
    ),
    "--config",
    str(p),
]
for p in files.values():
    cmd.extend(["--input", str(p)])
subprocess.run(cmd + ["--", "python", "scripts/focus_ft_primary.py"], cwd=root, check=True)
