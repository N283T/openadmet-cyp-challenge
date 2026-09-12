"""Summarize the completed report-derived round against the prior global pool."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


new, old = pl.read_csv(path(0)), pl.read_csv(path(1))
weights = json.loads(path(2).read_text())
catalog = {c["column"]: c for c in json.loads(path(3).read_text())}
summary = json.loads(path(4).read_text())
rows = []
for label, frame in [("previous_global", old), ("report_round_global", new)]:
    for r in frame.filter((pl.col("subset") == "all") & (pl.col("fold") == -1)).iter_rows(
        named=True
    ):
        rows.append(dict(round=label, **r))
comparison = pl.DataFrame(rows)
comparison.write_csv(out / "round-comparison.csv")
selected = []
for r in weights:
    for name, w in r["weights"].items():
        c = catalog[name]
        selected.append(
            dict(
                endpoint=r["endpoint"],
                method=r["method"],
                held_fold=r["held_fold"],
                candidate=name,
                weight=w,
                run=c["source"]["run"],
                label=c["label"],
                purpose=c["source"]["purpose"],
                new_in_round=c["source"]["run"] >= "20260912T092000",
            )
        )
s = pl.DataFrame(selected)
s.write_csv(out / "selected-models.csv")
paragraphs = [
    "# Report-derived experiments: final aggregate result",
    "",
    f"Expanded global selection `{cfg['selection_run']}`; independent audit `{cfg['audit_run']}`.",
    f"The pool expands from 107 to {summary['candidates']} candidate identities. "
    "All four research directions and a combined-descriptor check completed; "
    "model outputs are preserved.",
    "",
    "The comparison below uses leave-one-parent-fold-out weight fitting on the same "
    "retained OOF keys. Underlying models are reused: this remains a weight-stability "
    "check rather than newly nested model CV.",
    "",
    "| Target | Previous global | Expanded global | Difference |",
    "|---|---:|---:|---:|",
]
diffs = []
for ep in sorted(new["endpoint"].unique()):
    vals = []
    for frame in [old, new]:
        r = frame.filter(
            (pl.col("endpoint") == ep)
            & (pl.col("method") == "caruana_weight_crossfit")
            & (pl.col("subset") == "all")
            & (pl.col("fold") == -1)
        )
        vals.append(r["ST_RAE"].item())
    diffs.append(vals)
    paragraphs.append(
        f"| {ep.split('_')[0]} | {vals[0]:.6f} | {vals[1]:.6f} | {vals[1] - vals[0]:+.6f} |"
    )
before = sum(v[0] for v in diffs) / 4
after = sum(v[1] for v in diffs) / 4
paragraphs += [
    f"| Macro | {before:.6f} | {after:.6f} | {after - before:+.6f} |",
    "",
    "The macro movement is very small. CYP1A2 improves modestly; CYP2D6 barely moves. "
    "CYP2C9 and CYP3A4 regress slightly against the earlier global pool. No substantial "
    "focus-target breakthrough was found in this round.",
    "",
    "Full-OOF weight-fitting scores are listed separately because the same labels select "
    "the weights:",
]
for label, frame in [("Previous", old), ("Expanded", new)]:
    p = frame.filter(
        (pl.col("method") == "caruana_full_fit")
        & (pl.col("subset") == "all")
        & (pl.col("fold") == -1)
    )["ST_RAE"].mean()
    paragraphs.append(f"- {label}: {p:.6f}")
paragraphs += [
    "",
    "## New candidates selected in the full-OOF recipe",
    "",
    "This is the fitted recipe, not evidence of an independent gain for each selected "
    "component. Complete weights, including five separate weight-crossfit recipes per "
    "endpoint, are saved in selected-models.csv.",
]
for ep in sorted(new["endpoint"].unique()):
    paragraphs += [
        "",
        f"### {ep.split('_')[0]}",
        "",
        "| New candidate | Weight | Source |",
        "|---|---:|---|",
    ]
    h = s.filter(
        (pl.col("endpoint") == ep)
        & (pl.col("method") == "caruana_full_fit")
        & pl.col("new_in_round")
    ).sort("weight", descending=True)
    for r in h.iter_rows(named=True):
        description = (r["label"] + " — " + r["purpose"]).replace("|", "/")
        paragraphs.append(f"| `{r['run']}` | {r['weight']:.2%} | {description} |")
    if h.is_empty():
        paragraphs.append("| None | — | — |")
paragraphs += [
    "",
    "## Branch conclusions",
    "",
    "- TDI continuous potency and paired shifts produce small matched 1A2/2D6 LGBM gains "
    "but no macro gain in fixed additions.",
    "- Hill-linked TDI response heads underperform matched independent heads in the "
    "frozen-feature MLP; the MLP control itself is weaker than existing models. This is a "
    "negative result for this implementation, not a refutation of the mechanism.",
    "- Geometry features help CYP1A2 modestly, including the TabPFN comparison; CYP2D6 "
    "remains weak. Topology/basic-N and Gasteiger quantities are proxies, not "
    "pKa/microstate estimates.",
    "- Real external substrate probabilities have little downstream benefit. This branch "
    "tests substrate supervision; atom-level SoM remains untested.",
    "- Simultaneous geometry, substrate and Jazzy descriptors also fail to improve the "
    "prior manual ensemble at the prespecified fixed weights.",
    "",
    "All fitted models, descriptors and OOFs are retained. No full-data test refit, "
    "submission or production-pointer replacement is included in this round.",
]
(out / "report.md").write_text("\n".join(paragraphs) + "\n")
with psycopg.connect(**cfg["database"]) as conn:
    for r in new.iter_rows(named=True):
        conn.execute(
            "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
            (
                Path(os.environ["CYP_RUN_DIR"]).name,
                r["endpoint"],
                "report_round:" + r["subset"],
                str(r["fold"]),
                r["method"] + ":ST_RAE",
                r["ST_RAE"],
            ),
        )
print(before, after)
print(s.filter((pl.col("method") == "caruana_full_fit") & pl.col("new_in_round")))
