"""Independent fixed-five OOF reconstruction, official scoring and provenance package."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
fields = keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]


def read(name):
    return pl.read_parquet(inp / name)


base = read(cfg["final_oof"]).sort(keys)
assert base.height == 6384 and base.select(keys).n_unique() == 6384
blend_cfg = json.loads((inp / cfg["blend_config"]).read_text())
assert len(blend_cfg["components"]) == len(cfg["components"]) == 5
assert base["split_id"].unique().to_list() == [blend_cfg["split_run"]]
predictions, components = [], []
original_manifest = json.loads((inp / cfg["blend_manifest"]).read_text())
for i, component in enumerate(cfg["components"]):
    frame = read(component["oof"]).sort(keys)
    assert frame.height == 6384 and frame.select(keys).n_unique() == 6384
    for column in fields:
        assert frame[column].equals(base[column]), (i, column)
    column = blend_cfg["components"][i]["column"]
    assert blend_cfg["components"][i]["file"] == f"{i:03d}-external-oof.parquet"
    values = frame[column].to_numpy()
    assert np.isfinite(values).all()
    predictions.append(values)
    digest = hashlib.sha256((inp / component["oof"]).read_bytes()).hexdigest()
    assert digest == original_manifest["inputs"][i]["sha256"]
    result = json.loads((inp / component["result"]).read_text())
    components.append(
        component
        | dict(
            weight=0.2,
            prediction_sha256=digest,
            original_path=original_manifest["inputs"][i]["original"],
            run_result=result,
        )
    )
reconstructed = np.zeros(base.height)
for values in predictions:
    reconstructed += values
reconstructed /= 5
actual = base["prediction"].to_numpy()
assert np.isfinite(actual).all()
np.testing.assert_allclose(reconstructed, actual, rtol=0, atol=1e-12)
labels = read(cfg["measurements"]).filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
assignments = read(cfg["assignments"])
reference = (
    labels.join(assignments, on="compound_id", validate="m:1")
    .filter(~pl.col("is_anchor"))
    .sort(keys)
)
assert reference.height == 6384
for column in keys + ["conf_low", "conf_high", "is_neighbor"]:
    assert reference[column].equals(base[column]), column
assert reference["value"].equals(base["y_true"])
assert set(base["fold"]) == set(range(5))
assert (base["conf_low"] <= base["conf_high"]).all()
membership_checks = []
for filename in cfg["ft_membership"]:
    members = read(filename)
    fold = members["outer_fold"].unique().item()
    held = (members["fold"] == fold) & ~members["is_anchor"]
    assert (held == (members["role"] == "outer_holdout")).all()
    assert not members.filter(held & members["role"].is_in(["training", "inner_validation"])).height
    membership_checks.append(
        dict(file=filename, fold=fold, rows=members.height, excluded_nonanchors=int(held.sum()))
    )
nested_audit = json.loads((inp / cfg["nested_membership_audit"]).read_text())
assert len(nested_audit) == 10
assert len({tuple(r["excluded_folds"]) for r in nested_audit}) == 10
(out / "membership-verification.json").write_text(
    json.dumps(
        dict(actual_finetune=membership_checks, recorded_nested_audit=nested_audit), indent=2
    )
)
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
for (endpoint,), group in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        g = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = g if fold == -1 else g.filter(pl.col("fold") == fold)
            score = scorer.rae_soft_threshold_absolute_error(
                h["y_true"].to_numpy(),
                h["prediction"].to_numpy(),
                y_true_lower=h["conf_low"].to_numpy(),
                y_true_upper=h["conf_high"].to_numpy(),
            )
            rows.append(
                dict(endpoint=endpoint, subset=subset, fold=fold, n=h.height, ST_RAE=float(score))
            )
metrics = pl.DataFrame(rows).sort("subset", "endpoint", "fold")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("subset")
    .agg(pl.col("ST_RAE").mean(), pl.col("n").sum())
    .sort("subset")
)
prior = pl.read_csv(inp / cfg["prior_macro"]).filter(pl.col("method") == blend_cfg["method"])
for row in macro.to_dicts():
    expected = prior.filter(pl.col("subset") == row["subset"])["ST_RAE"].item()
    assert abs(expected - row["ST_RAE"]) < 1e-12
assert macro.filter(pl.col("subset") == "all")["ST_RAE"].item() < 0.6
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
base.with_columns(pl.Series("prediction", reconstructed)).write_parquet(
    out / "verified-oof.parquet"
)
hashes = [
    dict(snapshot=p.name, bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    for p in sorted(inp.iterdir())
    if p.is_file()
]
manifest = dict(
    config=cfg,
    components=components,
    input_hashes=hashes,
    bitwise_equal=bool(np.array_equal(actual, reconstructed)),
    max_absolute_difference=float(np.max(np.abs(actual - reconstructed))),
    rows=6384,
    official_macro=macro.to_dicts(),
    weight_fitting=False,
    upstream_checkpoint_role="Upstream 47-head pretraining, distinct from challenge 17-head FT",
    upstream_corpus_manifest="Unavailable; unknown composition alone is not evidence of leakage",
    scope="Independent OOF package only; no final all-training fit or test submission",
)
(out / "package-manifest.json").write_text(json.dumps(manifest, indent=2))
card = [
    "# Regression OOF result",
    "",
    "The fixed equal mean of five components achieves official macro ST-RAE below 0.6.",
    "",
    "| Subset | Rows | Macro ST-RAE |",
    "|---|---:|---:|",
]
for r in macro.to_dicts():
    card.append(f"| {r['subset']} | {r['n']} | {r['ST_RAE']:.12f} |")
card += ["", "## Components", "", "| Component | Run | Weight |", "|---|---|---:|"]
for c in components:
    card.append(f"| {c['name']} | {c['run_id']} | 0.2 |")
card += [
    "",
    "## Verification and provenance",
    "",
    "All 6,384 keyed compound/endpoint/fold rows, labels, confidence bounds and neighbor flags "
    "match the original measurements and retained parent assignments. Anchors are excluded from "
    "OOF evaluation. Original component hashes match the saved blend manifest. The independent "
    "equal mean reproduces predictions within 1e-12; official scores reproduce the prior audit.",
    "",
    "The user authorized public external pretrained weights. Jeremy upstream medium and CheMeleon "
    "checkpoints have 47 heads, distinct from his challenge 17-head fine-tune. "
    "Their original corpus manifest and named task schema are unavailable. This uncertainty "
    "does not itself establish leakage. Our downstream training uses the recorded outer-fold "
    "exclusions; detailed source manifests and audit references are included.",
    "",
    "These are development OOF results after multiple model comparisons, not an untouched test "
    "estimate or leaderboard result. The final ensemble uses fixed equal weights, with no OOF "
    "weight fit. This package does not contain a new all-training-data fit or submission.",
    "",
    "See package-manifest.json for every input hash and source reference; metrics.csv for "
    "target/fold/subset scores and verified-oof.parquet for keyed predictions.",
]
(out / "README.md").write_text("\n".join(card) + "\n")
print(
    json.dumps(
        dict(
            macro=macro.to_dicts(),
            bitwise_equal=manifest["bitwise_equal"],
            maximum_difference=manifest["max_absolute_difference"],
        ),
        indent=2,
    )
)
