"""Rebuild original UMAP log2fc auxiliaries with three folds excluded."""

import itertools
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = root / "artifacts/runs/shared/20260910T120234802703Z-540ba590"
base = json.loads((source / "config.json").read_text())
files = sorted((source / "inputs").iterdir())
assert len(files) == 5
directory = root / "configs/classification-stack-inner"
directory.mkdir(exist_ok=True)
registry = []
for triple in itertools.combinations(range(5), 3):
    cfg = dict(base, excluded_folds=list(triple), outer_fold=min(triple))
    config = directory / ("aux-" + "-".join(map(str, triple)) + ".json")
    config.write_text(json.dumps(cfg, indent=2))
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Historical classifier nested inner log2fc features: exclude UMAP triple {triple}",
        "--config",
        str(config),
    ]
    for path in files:
        command += ["--input", str(path)]
    child = subprocess.Popen(
        command + ["--", "python", "-m", "cyp_challenge.log2fc_pretrain"],
        cwd=root,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    run = None
    for line in child.stdout:
        print(line, end="", flush=True)
        if line.startswith("run_id="):
            run = line.strip().split("=", 1)[1]
    assert child.wait() == 0 and run
    registry.append(dict(excluded_folds=list(triple), run=run))
    (directory / "aux-registry.json").write_text(json.dumps(registry, indent=2))
