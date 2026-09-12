"""Verify expanded nested selection against original and report fixed OOF comparisons."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def table(key):
    return pl.read_parquet(inp / cfg[key])


def obj(key):
    return json.loads((inp / cfg[key]).read_text())


new, old = table("new_oof"), table("old_oof")
keys = ["compound_id", "endpoint", "fold"]
nm, om = table("new_membership"), table("old_membership")
assert nm.sort(keys + ["outer_fold"]).equals(om.sort(keys + ["outer_fold"]))
verification = []
for pair in cfg["inner_pairs"]:
    ni, oi = (pl.read_parquet(inp / pair[k]) for k in ["new", "old"])
    assert ni.width == oi.width + 6
    assert ni.select(oi.columns).equals(oi), pair
    verification.append({"new": pair["new"], "original_predictions_bitwise_equal": True})
preprocessor_checks = []
if cfg.get("classic"):
    classic = table("classic")
    cols = [c for c in classic.columns if c.startswith("mordred_")]
    for p in cfg["preprocessors"]:
        # These are our own recorded NumPy artifacts; compound IDs use object arrays.
        saved = np.load(inp / p["file"], allow_pickle=True)
        m = nm.filter(
            (pl.col("endpoint") == p["endpoint"]) & (pl.col("outer_fold") == p["outer_fold"])
        )
        m = m.filter(
            pl.col("role") == "inner_train"
            if p["role"] == "inner"
            else pl.col("role") != "outer_validation"
        )
        training_ids = saved["training_compound_ids"].tolist()
        assert sorted(training_ids) == sorted(m["compound_id"].to_list())
        x = (
            pl.DataFrame({"compound_id": training_ids})
            .join(classic, on="compound_id", validate="1:1", maintain_order="left")
            .select(cols)
            .to_numpy()
            .astype(np.float32)
        )
        x[~np.isfinite(x)] = np.nan
        keep = np.flatnonzero(np.isfinite(x).any(axis=0))
        keep = keep[np.nanmin(x[:, keep], axis=0) != np.nanmax(x[:, keep], axis=0)]
        assert np.array_equal(keep, saved["descriptor_indices"])
        assert np.array_equal(np.nanmedian(x[:, keep], axis=0), saved["medians"])
        preprocessor_checks.append(
            dict(
                **p,
                verified_training_ids=True,
                verified_filter_and_medians=True,
                retained_descriptors=len(keep),
            )
        )
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
baseline = table("baseline").sort(keys)
new, old, jeremy = (t.sort(keys) for t in [new, old, table("jeremy_oof")])
ft = table("jeremy_ft_oof").sort(keys)
for t in [new, old, jeremy, ft]:
    assert t.height == 6384 and t.select(keys).equals(baseline.select(keys))
    for c in ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert t[c].equals(baseline[c]), c
base = baseline.with_columns(
    new["prediction"].alias("expanded"),
    old["prediction"].alias("original"),
    jeremy["prediction"].alias("frozen_jeremy_pfn"),
    ft["prediction"].alias("jeremy_ft_logfc_lgbm"),
)
base = base.with_columns(
    pl.mean_horizontal("expanded", "frozen_jeremy_pfn").alias("expanded_jeremy_fixed_pair"),
    pl.mean_horizontal("original", "frozen_jeremy_pfn").alias("original_jeremy_fixed_pair"),
    pl.mean_horizontal("expanded", "frozen_jeremy_pfn", "jeremy_ft_logfc_lgbm").alias(
        "expanded_fixed_three"
    ),
    pl.mean_horizontal("original", "frozen_jeremy_pfn", "jeremy_ft_logfc_lgbm").alias(
        "original_fixed_three"
    ),
)
methods = [
    "expanded",
    "original",
    "frozen_jeremy_pfn",
    "expanded_jeremy_fixed_pair",
    "original_jeremy_fixed_pair",
    "mean_four_pubchem600",
    "expanded_fixed_three",
    "original_fixed_three",
]
if cfg.get("previous_expanded_oof"):
    prev = table("previous_expanded_oof").sort(keys)
    assert prev.select(keys).equals(base.select(keys))
    base = base.with_columns(prev["prediction"].alias("previous_expanded"))
    methods.append("previous_expanded")
metrics = []
for (ep,), group in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        h = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for method in methods:
            metrics.append(
                dict(
                    endpoint=ep,
                    subset=subset,
                    method=method,
                    ST_RAE=float(
                        scorer.rae_soft_threshold_absolute_error(
                            h["y_true"].to_numpy(),
                            h[method].to_numpy(),
                            y_true_lower=h["conf_low"].to_numpy(),
                            y_true_upper=h["conf_high"].to_numpy(),
                        )
                    ),
                )
            )
metrics = pl.DataFrame(metrics)
macro = metrics.group_by("method", "subset").agg(pl.col("ST_RAE").mean()).sort("subset", "ST_RAE")
metrics.write_csv(out / "comparison.csv")
macro.write_csv(out / "macro.csv")
base.write_parquet(out / "aligned-oof.parquet")
weights = []
for fit in obj("new_summary")["selections"]:
    for c in fit["selected"]:
        weights.append(
            dict(
                endpoint=fit["endpoint"],
                fold=fit["fold"],
                view=c["view"],
                loss=c["params"]["objective"],
                alpha=c["params"].get("alpha"),
                weight=c["weight"],
            )
        )
w = pl.DataFrame(weights, infer_schema_length=None)
w.write_csv(out / "selected-weights.csv")
family = (
    w.group_by("view")
    .agg(
        pl.col("weight").sum() / 20, pl.struct("endpoint", "fold").n_unique().alias("selected_fits")
    )
    .sort("weight", descending=True)
)
family.write_csv(out / "family-frequency.csv")
w.group_by("loss", "alpha").agg(
    pl.col("weight").sum() / 20, pl.struct("endpoint", "fold").n_unique().alias("selected_fits")
).write_csv(out / "loss-frequency.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            memberships_identical=True,
            original_candidate_checks=verification,
            preprocessor_checks=preprocessor_checks,
            matched_oof=6384,
            weight_optimization_on_oof=False,
            upstream_pretrained_diagnostic=True,
        ),
        indent=2,
    )
)
lines = [
    "# " + cfg.get("title", "Jeremy-expanded nested Caruana diagnostic"),
    "",
    f"Fit `{cfg['fit_run']}` uses {cfg['original_candidate_count'] + 6} candidate "
    "LightGBMs with L1, quantile 0.35 and quantile 0.65 losses. "
    f"The previous {cfg['original_candidate_count']} candidate "
    "inner predictions and all membership assignments are bitwise identical to the earlier run. "
    "Each of 20 endpoint/fold fits selects 50 replacement steps using its 20% inner group holdout, "
    "retains the best inner iteration and refits the selected models on all outer training labels.",
    "",
    "The frozen upstream encoder is a separate 47-head pretraining checkpoint. "
    "Its training molecule manifest is unavailable; these scores are diagnostic "
    "and are not adopted as the formal goal result.",
    "",
    "All 6,384 OOF identities, endpoints, folds, labels and bounds match the fixed baseline. "
    "The two fixed pairs use equal 0.5 weights. The fixed three-way blends add "
    "Jeremy fine-tuned embedding + log2fc LGBM with equal thirds. "
    "No OOF weight optimization is performed.",
    "For FP+Mordred variants, all 40 saved preprocessing artifacts were independently "
    "checked: exact inner/outer training IDs, finite/nonconstant descriptor filters and "
    "training medians. The original score comparator below remains the 24-candidate ensemble."
    if preprocessor_checks
    else "",
    "",
    "| Method | All ST-RAE | Neighbors ST-RAE |",
    "|---|---:|---:|",
]
for method in methods:
    v = {r["subset"]: r["ST_RAE"] for r in macro.filter(pl.col("method") == method).to_dicts()}
    lines.append(f"| {method} | {v['all']:.6f} | {v['neighbors']:.6f} |")
lines += [
    "",
    "## Endpoint results",
    "",
    "| CYP | Expanded | Original | Expanded + Jeremy PFN | Original + Jeremy PFN |",
    "|---|---:|---:|---:|---:|",
]
for ep in sorted(metrics["endpoint"].unique()):
    v = {
        r["method"]: r["ST_RAE"]
        for r in metrics.filter((pl.col("endpoint") == ep) & (pl.col("subset") == "all")).to_dicts()
    }
    lines.append(
        f"| {ep.split('_')[0]} | {v['expanded']:.6f} | {v['original']:.6f} | "
        f"{v['expanded_jeremy_fixed_pair']:.6f} | {v['original_jeremy_fixed_pair']:.6f} |"
    )
lines += [
    "",
    "## Inner-selected families",
    "",
    "Weights average across all 20 fits, including zero weights.",
    "",
    "| View | Mean weight | Selected fits |",
    "|---|---:|---:|",
]
for r in family.to_dicts():
    lines.append(f"| {r['view']} | {r['weight']:.4f} | {r['selected_fits']} |")
(out / "report.md").write_text("\n".join(lines) + "\n")
print(macro, flush=True)
print(family, flush=True)
