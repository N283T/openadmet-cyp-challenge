"""Launch fixed, equal-weight TDI ensembles from audited inner and outer fits."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(run, file):
    p = root / "artifacts/runs/classification" / run / "outputs" / file
    assert p.is_file(), p
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


lg = "20260912T115220209511Z-76facc71"
pf = "20260912T115242271820Z-51023dac"
graph = "20260912T120135164747Z-1da37202"
sources = [
    dict(
        prefix="lgbm",
        oof=add(lg, "oof.parquet"),
        inner=add("20260912T120959802121Z-8aeb9a25", "inner-predictions.parquet"),
        test=add(lg, "fold-test-probabilities.parquet"),
    ),
    dict(
        prefix="pfn",
        oof=add(pf, "oof-probabilities.parquet"),
        inner=add("20260912T115749418593Z-397ddee7", "inner-validation-probabilities.parquet"),
        test=add(pf, "test-fold-probabilities.parquet"),
    ),
    dict(
        prefix="graph",
        oof=add(graph, "classification-oof.parquet"),
        inner=add(graph, "inner-probabilities.parquet"),
        test=add(graph, "fold-test-probabilities.parquet"),
    ),
]
views = [
    "morgan",
    "mordred",
    "fp_mordred",
    "external_class",
    "external_mixed",
    "external_class_mordred",
]
plain = [f"lgbm::{v}_plain" for v in views]
balanced = [f"lgbm::{v}_balanced" for v in views]
pfn = [f"pfn::{v}_pca256" for v in ["external_class", "external_mixed", "chemeleon"]]
graphs = ["graph::tdi_only", "graph::tdi_aux"]
cfg = dict(
    sources=sources,
    ensembles={
        "classical_plain_mean": plain[:3],
        "lgbm_plain_mean": plain,
        "lgbm_balanced_mean": balanced,
        "pfn_mean": pfn,
        "graph_mean": graphs,
        "classical_pfn_mean": plain[:3] + pfn,
        "all_plain_pfn_graph_mean": plain + pfn + graphs,
    },
)
config = root / "configs/classification-round/fixed-ensemble-v2.json"
assert not config.exists()
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Compare seven fixed probability means using isolated inner thresholds",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_fixed_ensemble.py"], cwd=root, check=True
)
