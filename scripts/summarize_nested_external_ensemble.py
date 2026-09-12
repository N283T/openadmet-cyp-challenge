"""Summarize fixed audit metrics and training-only Caruana selection frequencies."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
a = pl.read_csv(inp / "000-metrics.csv")
b = pl.read_csv(inp / "001-metrics.csv")
s = json.loads((inp / "002-summary.json").read_text())
methods = [
    "nested_caruana_external_lgbm",
    "nested_caruana_external_lgbm__four",
    "mean_four_pubchem600",
]
best = "external3_class_mt_v1_ft_tabpfn__four"
assert best in set(b["method"])
comp = pl.concat([a.filter(pl.col("method").is_in(methods)), b.filter(pl.col("method") == best)])
comp = comp.filter(pl.col("fold") == -1).sort("subset", "endpoint", "method")
comp.write_csv(out / "comparison.csv")
macro = comp.group_by("method", "subset").agg(pl.col("ST_RAE").mean()).sort("subset", "ST_RAE")
macro.write_csv(out / "macro.csv")
weights = []
for record in s["selections"]:
    for c in record["selected"]:
        weights.append(
            dict(
                endpoint=record["endpoint"],
                fold=record["fold"],
                view=c["view"],
                objective=c["params"]["objective"],
                alpha=c["params"].get("alpha"),
                weight=c["weight"],
            )
        )
w = pl.DataFrame(weights, infer_schema_length=None)
w.write_csv(out / "selected-weights.csv")
family = (
    w.group_by("view")
    .agg(
        pl.col("weight").sum().alias("total_weight"),
        pl.struct("endpoint", "fold").n_unique().alias("selected_fits"),
    )
    .with_columns((pl.col("total_weight") / 20).alias("mean_weight"))
    .sort("mean_weight", descending=True)
)
family.write_csv(out / "family-frequency.csv")
loss = (
    w.group_by("objective", "alpha")
    .agg(
        pl.col("weight").sum().alias("total_weight"),
        pl.struct("endpoint", "fold").n_unique().alias("selected_fits"),
    )
    .with_columns((pl.col("total_weight") / 20).alias("mean_weight"))
    .sort("mean_weight", descending=True)
)
loss.write_csv(out / "loss-frequency.csv")
lines = [
    "# Nested external Caruana ensemble",
    "",
    "Run `20260912T054406239560Z-96fce32b` fitted 24 candidate LightGBMs per endpoint and "
    "outer fold: eight feature views crossed with L1 and quantile losses at 0.35 and 0.65. "
    "Each used 300 trees, 15 leaves, min_child_samples 20, learning rate 0.05 and L2=1. "
    "Three CPU threads; no GPU.",
    "",
    "Weights were selected on a 20% inner group holdout of outer-training compounds, "
    "retaining all anchor-containing groups in inner training. Caruana ran 50 replacement "
    "steps and retained the iteration with lowest inner official ST-RAE. Selected candidates "
    "were refitted on all outer-training observations. Outer-fold labels were not used "
    "to fit weights or choose iterations.",
    "",
    "All 6,384 OOF predictions matched the baseline identities, endpoints, folds and bounds. "
    "Audit `20260912T055242888594Z-56ff81b1` evaluates the nested candidate plus the "
    "already prescribed equal blends; no OOF weight optimization was performed.",
    "",
    "## Results",
    "",
    "| Method | All OOF ST-RAE | Neighbor ST-RAE |",
    "|---|---:|---:|",
]
for method in methods + [best]:
    v = {r["subset"]: r["ST_RAE"] for r in macro.filter(pl.col("method") == method).to_dicts()}
    lines.append(f"| {method} | {v['all']:.6f} | {v['neighbors']:.6f} |")
lines += [
    "",
    "Four means an equal blend of the candidate, FP+Mordred LGBM, predicted-log2fc TabPFN "
    "and frozen Admetica-embedding TabPFN. The original baseline is FP+Mordred, "
    "predicted-log2fc, Admetica and PubChem600. The FT reference comes from "
    "audit `20260912T054319015342Z-f96a88a3`.",
    "",
    "The nested ensemble improves the original baseline, especially on the neighbor subset, "
    "but its best prescribed blend does not beat the current FT-embedding TabPFN "
    "four-model blend overall. It is retained as a diverse candidate, "
    "not adopted as a new global best.",
    "",
    "## Endpoint comparison: all OOF",
    "",
    "| CYP | Original baseline | Nested Caruana | Caruana four | FT TabPFN four |",
    "|---|---:|---:|---:|---:|",
]
for ep in sorted(comp["endpoint"].unique()):
    vals = {
        r["method"]: r["ST_RAE"]
        for r in comp.filter((pl.col("endpoint") == ep) & (pl.col("subset") == "all")).to_dicts()
    }
    lines.append(
        f"| {ep.split('_')[0]} | {vals[methods[2]]:.6f} | {vals[methods[0]]:.6f} | "
        f"{vals[methods[1]]:.6f} | {vals[best]:.6f} |"
    )
lines += [
    "",
    "## Inner-selected feature families",
    "",
    "Mean weights average across all 20 endpoint/fold fits, including zero weights. "
    "Selected fits counts a family once per fit even when multiple losses were selected. "
    "These are training-selection frequencies, not a causal feature-importance measure.",
    "",
    "| View | Mean weight | Selected fits /20 |",
    "|---|---:|---:|",
]
for r in family.to_dicts():
    lines.append(f"| {r['view']} | {r['mean_weight']:.4f} | {r['selected_fits']} |")
lines += ["", "| Loss | Quantile | Mean weight | Selected fits /20 |", "|---|---:|---:|---:|"]
for r in loss.to_dicts():
    lines.append(
        f"| {r['objective']} | {r['alpha']} | {r['mean_weight']:.4f} | {r['selected_fits']} |"
    )
(out / "report.md").write_text("\n".join(lines) + "\n")
print(macro, flush=True)
print(comp.filter(pl.col("subset") == "all").select("endpoint", "method", "ST_RAE"), flush=True)
