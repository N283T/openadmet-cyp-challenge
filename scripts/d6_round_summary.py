"""Summarize fixed-prior versus expanded D6 pool, preserving fitted/CV distinction."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


new, old = pl.read_csv(path(0)), pl.read_csv(path(1))
weights = json.loads(path(2).read_text())
catalog = {r["column"]: r for r in json.loads(path(3).read_text())}
summary = json.loads(path(4).read_text())
previous, pool, current = [pl.read_parquet(path(i)) for i in range(5, 8)]
spec = importlib.util.spec_from_file_location("official", path(8))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
comparisons = []
for label, frame in [("previous", old), ("expanded", new)]:
    for r in frame.filter((pl.col("subset") == "all") & (pl.col("fold") == -1)).iter_rows(
        named=True
    ):
        comparisons.append(dict(round=label, **r))
pl.DataFrame(comparisons).write_csv(out / "global-comparison.csv")
selected = []
for w in weights:
    for name, value in w["weights"].items():
        c = catalog[name]
        selected.append(
            dict(
                endpoint=w["endpoint"],
                method=w["method"],
                held_fold=w["held_fold"],
                candidate=name,
                weight=value,
                run=c["source"]["run"],
                label=c["label"],
                new_in_round=(
                    c["source"]["run"].startswith("2026")
                    and c["source"]["run"] >= "20260912T100000"
                ),
            )
        )
pl.DataFrame(selected).write_csv(out / "selected-models.csv")
ep = "CYP2D6_pIC50_direct_inhibition"
h = pool.filter(pl.col("endpoint") == ep).sort(keys)
b = previous.filter(
    (pl.col("endpoint") == ep) & (pl.col("method") == "caruana_weight_crossfit")
).sort(keys)
for c in keys + ["y_true", "conf_low", "conf_high"]:
    assert h[c].equals(b[c])


def score(t, p):
    return float(
        scorer.rae_soft_threshold_absolute_error(
            t["y_true"].to_numpy(),
            p,
            y_true_lower=t["conf_low"].to_numpy(),
            y_true_upper=t["conf_high"].to_numpy(),
        )
    )


standalone = []
for name, c in catalog.items():
    if (
        not c["source"]["run"].startswith("2026")
        or c["source"]["run"] < "20260912T100000"
        or ep not in c["endpoints"]
    ):
        continue
    p = h[name].to_numpy()
    for weight in [0.25, 0.5, 1.0]:
        pred = (1 - weight) * b["prediction"].to_numpy() + weight * p
        standalone.append(
            dict(
                run=c["source"]["run"],
                label=c["label"],
                weight=weight,
                ST_RAE=score(h, pred),
                Spearman=float(spearmanr(h["y_true"], pred).statistic),
            )
        )
pl.DataFrame(standalone).write_csv(out / "new-candidate-comparison.csv")
lines = [
    "# CYP2D6 focused round: aggregate result",
    "",
    f"Selection `{cfg['selection_run']}`; audit `{cfg['audit_run']}`.",
    f"Expanded pool: {summary['candidates']} candidate identities.",
    "",
    "The comparison holds each parent fold out when fitting ensemble weights.",
    "Base OOF models are reused, so this checks weight stability, not fully nested model CV.",
    "",
    "| Endpoint | Previous | Expanded | Difference |",
    "|---|---:|---:|---:|",
]
for endpoint in sorted(previous["endpoint"].unique()):
    left = previous.filter(
        (pl.col("endpoint") == endpoint) & (pl.col("method") == "caruana_weight_crossfit")
    ).sort(keys)
    right = current.filter(
        (pl.col("endpoint") == endpoint) & (pl.col("method") == "caruana_weight_crossfit")
    ).sort(keys)
    a, z = score(left, left["prediction"].to_numpy()), score(right, right["prediction"].to_numpy())
    lines.append(f"| {endpoint} | {a:.6f} | {z:.6f} | {z - a:+.6f} |")
    if endpoint != ep:
        assert (
            np.max(np.abs(left["prediction"].to_numpy() - right["prediction"].to_numpy())) < 1e-12
        ), endpoint
lines += [
    "",
    "Other endpoint predictions match the frozen prior pool within 1e-12.",
    "",
    "## New components in full-OOF fitted CYP2D6 recipe",
    "",
    "These weights use the full OOF labels and are not an independent validation result.",
    "",
    "| Run | Method | Weight |",
    "|---|---|---:|",
]
for s in selected:
    if s["endpoint"] == ep and s["method"] == "caruana_full_fit" and s["new_in_round"]:
        lines.append(f"| `{s['run']}` | {s['label']} | {s['weight']:.2%} |")
(out / "report.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines), flush=True)
