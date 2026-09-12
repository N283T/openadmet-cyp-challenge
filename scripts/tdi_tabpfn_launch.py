"""Snapshot frozen external feature inputs for the bounded TDI TabPFN comparison."""

import json
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
prep = runs / "classification/20260912T115018743045Z-269e598d/outputs"
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


cfg = dict(
    labels=add(prep / "classification-labels.parquet"),
    assignments=add(prep / "umap-assignments.parquet"),
    test_mapping=add(runs / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"),
    weights=add(
        runs
        / "shared/20260910T121624900089Z-271684c2/inputs"
        / "004-tabpfn-v2.6-classifier-v2.6_default.ckpt"
    ),
    expected_counts={"CYP2D6_is_TDI": 1497, "CYP3A4_is_TDI": 3584},
    split_id="original_umap_269e598d",
    views=[],
)
for name, run, pattern, prefix in [
    (
        "external_class_pca256",
        "20260912T051543210441Z-a3ea2ec3",
        "external-embeddings.parquet",
        "emb_",
    ),
    (
        "external_mixed_pca256",
        "20260912T052853359792Z-2c0a0953",
        "external-embeddings.parquet",
        "emb_",
    ),
    ("chemeleon_pca256", "20260910T012630927160Z-e9576045", "chemeleon-*.parquet", "chemeleon_"),
]:
    feature_paths = sorted((runs / "shared" / run / "outputs").glob(pattern))
    assert feature_paths
    cfg["views"].append(
        dict(name=name, source_run=run, prefix=prefix, files=[add(path) for path in feature_paths])
    )
inner = "--inner" in sys.argv
if inner:
    outer_run = "20260912T115242271820Z-51023dac"
    cfg["outer_run"] = outer_run
    cfg["outer_oof"] = add(
        runs / "classification" / outer_run / "outputs/oof-probabilities.parquet"
    )
    cfg["outer_test"] = add(
        runs / "classification" / outer_run / "outputs/test-fold-probabilities.parquet"
    )
config = root / "configs/classification-tabpfn" / ("inner-v1.json" if inner else "frozen-v1.json")
config.parent.mkdir(exist_ok=True)
config.write_text(json.dumps(cfg, indent=2))
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Thirty inner UMAP threshold fits" if inner else "Thirty fixed UMAP TDI TabPFN2.6 fits",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
script = "scripts/tdi_tabpfn_inner.py" if inner else "scripts/tdi_tabpfn_frozen.py"
subprocess.run(command + ["--", "python", script], cwd=root, check=True)
