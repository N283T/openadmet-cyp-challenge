"""Save a readable global-selection report and register audited metrics."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


catalog = {x["column"]: x for x in json.loads(path(0).read_text())}
weights = json.loads(path(1).read_text())
metrics, macro = pl.read_csv(path(2)), pl.read_csv(path(3))
summary = json.loads(path(4).read_text())
selected = []
for row in weights:
    for col, w in row["weights"].items():
        c = catalog[col]
        selected.append(
            dict(
                endpoint=row["endpoint"],
                method=row["method"],
                held_fold=row["held_fold"],
                candidate=col,
                weight=w,
                run=c["source"]["run"],
                script=c["source"]["script"],
                label=c["label"],
                purpose=c["source"]["purpose"],
            )
        )
selection = pl.DataFrame(selected)
selection.write_csv(out / "selected-models.csv")
text = [
    "# Global Caruana across saved regression model families",
    "",
    "The user requested actual selection across all available compatible model families, "
    "replacing the prior manual top-level average.",
    "",
    f"Selection run: `{cfg['selection_run']}`. Independent audit: `{cfg['audit_run']}`.",
    "",
    f"The frozen inventory contains {summary['sources']} completed OOF artifacts. "
    "Exact endpoint/fold/label alignment and duplicate-vector removal leave "
    f"{summary['candidates']} candidate identities "
    f"({min(summary['endpoints'].values())}–{max(summary['endpoints'].values())} "
    "eligible per endpoint).",
    "Candidates include LightGBM, TabPFN, TabICL, Jazzy additions, FP/Mordred, external "
    "representations, neural direct predictions, and previously nested Caruana components. "
    "Per-fold early FP/Mordred and separate embedding/log2fc readouts are assembled explicitly. "
    "The pool also retains prior fixed-ensemble predictions exposed by independent final-package "
    "verification runs (m0079 and m0080), alongside individual model candidates. Thus the pool "
    "contains both individual models and aggregate predictions; it is not exclusively atomic "
    "models. Calibration and other comparison scripts are excluded by the recorded inventory "
    "filter, but this is not a claim that every earlier manual blend was excluded. Exact "
    "duplicate vectors are collapsed. All source exclusions and duplicate aliases are saved.",
    "",
    "For each endpoint, greedy selection with replacement runs 50 steps and retains the best "
    "training-loss prefix. Interval MAE is minimized; the ST-RAE denominator is fixed within "
    "the endpoint, so this optimizes the official ST-RAE objective. No candidate score "
    "prefilter, new hyperparameter search, or clipping is used.",
    "",
    "Two outputs are distinct: full-OOF fitted weights (development fit), and "
    "leave-one-parent-fold-out weight fitting. The latter excludes evaluation-fold labels "
    "from weight fitting only. Base models were not retrained in newly nested splits, and "
    "some selector-training predictions come from models trained using that evaluation fold. "
    "Therefore it is a weight-crossfit sensitivity check, not an independent nested-CV estimate.",
    "",
    "## Scores",
    "",
    "| Endpoint | Previous candidate | Weight cross-fit | Full-OOF weight fit |",
    "|---|---:|---:|---:|",
]
for ep in sorted(metrics["endpoint"].unique()):
    m = metrics.filter(
        (pl.col("endpoint") == ep) & (pl.col("subset") == "all") & (pl.col("fold") == -1)
    )
    vals = {r["method"]: r["ST_RAE"] for r in m.iter_rows(named=True)}
    text.append(
        f"| {ep.split('_')[0]} | {vals['latest']:.6f} | "
        f"{vals['caruana_weight_crossfit']:.6f} | {vals['caruana_full_fit']:.6f} |"
    )
vals = {
    r["method"]: r["ST_RAE"] for r in macro.filter(pl.col("subset") == "all").iter_rows(named=True)
}
text += [
    f"| Macro | {vals['latest']:.6f} | {vals['caruana_weight_crossfit']:.6f} | "
    f"{vals['caruana_full_fit']:.6f} |",
    "",
    "Global weight selection improves the macro point estimate chiefly through 2C9 and 3A4. "
    "The weight-crossfit results slightly worsen both focus endpoints, 1A2 and 2D6. This is "
    "not evidence that their underlying model problem has been solved. Neighbor macro is "
    "also essentially unchanged/slightly worse (0.675304 to 0.676179). MAE and Spearman "
    "are preserved in the audit and do not improve uniformly.",
    "",
    "Retain the global weights, aligned pool, selection histories, correlations, and keyed "
    "predictions as a new ensemble experiment. The lower full-fit score must not be "
    "presented as an untouched CV gain. No production pointer, full-data test refit, "
    "or submission was changed.",
    "",
    "## Full-OOF fitted composition",
    "",
    "Weights below are the full-OOF selection recipe, not a new independently validated "
    "deployment. Exact cross-fitted weights are in `selected-models.csv` and `weights.json`.",
]
for ep in sorted(metrics["endpoint"].unique()):
    text += ["", f"### {ep.split('_')[0]}", "", "| Candidate | Weight | Source |", "|---|---:|---|"]
    h = selection.filter(
        (pl.col("endpoint") == ep) & (pl.col("method") == "caruana_full_fit")
    ).sort("weight", descending=True)
    for r in h.iter_rows(named=True):
        source = (r["label"] + " — " + r["purpose"]).strip(" —").replace("|", "/")
        text.append(f"| `{r['candidate']}` / `{r['run']}` | {r['weight']:.2%} | {source} |")
text += [
    "",
    "## Verification",
    "",
    "The independent audit reconstructs all predictions from nonnegative, sum-one weights; "
    "verifies 6,384 exact keys and original confidence bounds; checks all 24 weight sets "
    "and selection-fold masks; and recomputes official pooled, fold, and neighbor metrics. "
    "Maximum reconstruction error is below 2e-15. The previous 0.5869306288697727 score "
    "is reproduced.",
    "",
    "Saved artifacts: `catalog.json`, `excluded.json`, `aliases.json`, `aligned-pool.parquet`, "
    "`global-oof.parquet`, `weights.json`, `selection-history.json`, `selected-models.csv`, "
    "and independent `metrics.csv`/`macro.csv`.",
]
(out / "report.md").write_text("\n".join(text) + "\n")
with psycopg.connect(**cfg["database"]) as conn:
    for r in metrics.iter_rows(named=True):
        conn.execute(
            "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
            (
                os.environ["CYP_RUN_DIR"].split("/")[-1],
                r["endpoint"],
                "global_caruana:" + r["subset"],
                str(r["fold"]),
                r["method"] + ":ST_RAE",
                r["ST_RAE"],
            ),
        )
print(
    selection.filter(pl.col("method") == "caruana_full_fit").select(
        "endpoint", "candidate", "weight", "run", "label"
    )
)
