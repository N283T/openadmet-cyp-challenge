"""Paired outer-evaluation comparisons and deployed ensemble complexity."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import matthews_corrcoef

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["oof"])
ref = pl.read_parquet(inp / cfg["reference"])
deploy = pl.read_parquet(inp / cfg["deployment"])
selections = json.loads((inp / cfg["selections"]).read_text())
metrics = pl.read_csv(inp / cfg["metrics"])
labels = pl.read_parquet(inp / cfg["labels"]).select("compound_id", "endpoint", "group")
rng = np.random.default_rng(1729)
rows, weights = [], []
for (ep, method), t in oof.group_by("endpoint", "method"):
    t = t.sort("compound_id").join(
        labels, on=["compound_id", "endpoint"], validate="1:1", maintain_order="left"
    )
    rname = "all_plain_pfn_graph_mean" if ep == "CYP2D6_is_TDI" else "lgbm_balanced_mean"
    r = ref.filter((pl.col("endpoint") == ep) & (pl.col("method") == rname)).sort("compound_id")
    assert t["compound_id"].equals(r["compound_id"])
    assert np.array_equal(t["y_true"], r["y_true"])
    y = t["y_true"].to_numpy()
    p = (t["prediction"] >= t["threshold"]).to_numpy()
    rp = (r["prediction"] >= r["threshold"]).to_numpy()
    group = t["group"].to_numpy()
    strata = [np.flatnonzero(group == g) for g in np.unique(group)]
    boot = []
    for _ in range(500):
        ix = np.concatenate([strata[i] for i in rng.integers(0, len(strata), len(strata))])
        if len(np.unique(y[ix])) < 2:
            continue
        boot.append(matthews_corrcoef(y[ix], p[ix]) - matthews_corrcoef(y[ix], rp[ix]))
    chosen = [s for s in selections if s["endpoint"] == ep and s["recipe"] == method]
    final = next(s for s in chosen if s["fold"] == -1)
    active = {k: v for k, v in final["weights"].items() if v > 0}
    for s in chosen:
        for name, weight in s["weights"].items():
            if weight > 0:
                weights.append(
                    dict(
                        endpoint=ep,
                        recipe=method,
                        fold=s["fold"],
                        model=name,
                        weight=weight,
                        threshold=s["threshold"],
                    )
                )
    blind = deploy.filter((pl.col("endpoint") == ep) & (pl.col("method") == method))
    rows.append(
        dict(
            endpoint=ep,
            method=method,
            reference=rname,
            MCC=float(matthews_corrcoef(y, p)),
            reference_MCC=float(matthews_corrcoef(y, rp)),
            difference=float(matthews_corrcoef(y, p) - matthews_corrcoef(y, rp)),
            conditional_group_lo=float(np.quantile(boot, 0.025)),
            conditional_group_hi=float(np.quantile(boot, 0.975)),
            final_base_types=len(active),
            final_readout_instances=5 * len(active),
            final_cut=final["threshold"],
            test_positives=int(blind["binary"].sum()),
            outer_types_min=min(
                sum(v > 0 for v in s["weights"].values()) for s in chosen if s["fold"] >= 0
            ),
            outer_types_max=max(
                sum(v > 0 for v in s["weights"].values()) for s in chosen if s["fold"] >= 0
            ),
        )
    )
pl.DataFrame(rows).write_csv(out / "paired-comparison.csv")
pl.DataFrame(weights).write_csv(out / "active-weights.csv")
metrics.filter(pl.col("fold") == -1).write_csv(out / "outer-metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            uncertainty="500 paired UMAP-group draws; conditional on fitted/screened models",
            complexity="Fivefold readouts; historical auxiliary models are additional dependencies",
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows).write_csv())
