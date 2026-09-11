"""Paired untouched-fold evaluation of auxiliary HPO and downstream regression."""

import importlib.util
import json
import os
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl
import psycopg
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
frames = {}
for method in ["original", "tuned"]:
    frames[method] = (
        pl.read_parquet(inp / cfg[method])
        .filter((pl.col("mode") == "parent_retained") & (pl.col("features") == "full"))
        .sort("endpoint", "compound_id", "fold")
    )
keys = ["endpoint", "compound_id", "fold"]
assert frames["original"].select(keys).equals(frames["tuned"].select(keys))
for col in ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert frames["original"][col].equals(frames["tuned"][col])
rows = []
for method, frame in frames.items():
    for (endpoint,), group in frame.group_by("endpoint"):
        for subset in ["all", "neighbors"]:
            h = group if subset == "all" else group.filter(pl.col("is_neighbor"))
            y, p, lo, hi = [
                h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]
            ]
            rows.append(
                {
                    "method": method,
                    "endpoint": endpoint,
                    "subset": subset,
                    "n": h.height,
                    "ST_RAE": float(
                        official.rae_soft_threshold_absolute_error(
                            y, p, y_true_lower=lo, y_true_upper=hi
                        )
                    ),
                    "MAE": float(np.abs(y - p).mean()),
                    "Spearman": float(spearmanr(y, p).statistic),
                }
            )
metrics = pl.DataFrame(rows).sort("subset", "endpoint", "method")
metrics.write_csv(out / "readout-metrics.csv")
with psycopg.connect(**cfg["database"]) as conn:
    for r in metrics.iter_rows(named=True):
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "rustuna_holdout:" + r["subset"],
                    "0",
                    r["method"] + ":" + metric,
                    r[metric],
                ),
            )
macro = (
    metrics.group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "method")
)
macro.write_csv(out / "readout-macro.csv")
assignments = pl.read_parquet(inp / cfg["assignments"]).filter(pl.col("fold") == 0)
targets = pl.read_parquet(inp / cfg["targets"])
aux_rows, aligned = [], []
for method in ["original", "tuned"]:
    ps = []
    for paths in cfg["aux"][method]:
        f = pl.concat([pl.read_parquet(inp / p) for p in paths])
        ps.append(f.unpivot(index="compound_id", variable_name="head", value_name="pred_log2fc"))
    mean = (
        pl.concat(ps)
        .group_by("compound_id", "head")
        .agg(pl.col("pred_log2fc").mean(), pl.len().alias("models"))
    )
    assert (mean["models"] == 4).all()
    f = mean.join(assignments.select("compound_id", "fold"), on="compound_id", validate="m:1").join(
        targets, on=["compound_id", "head"], validate="1:1"
    )
    for (head,), h in f.group_by("head"):
        aux_rows.append(
            {
                "method": method,
                "head": head,
                "scope": "all_primary_fold0",
                "n": h.height,
                "RMSE": float(np.sqrt(((h["pred_log2fc"] - h["value"]) ** 2).mean())),
                "Spearman": float(spearmanr(h["pred_log2fc"], h["value"]).statistic),
            }
        )
    pred = mean.with_columns(
        (pl.col("head").str.split("@").list.first() + "_pIC50_direct_inhibition").alias("endpoint")
    )
    h = frames[method].join(
        pred.select("compound_id", "endpoint", "pred_log2fc"),
        on=["compound_id", "endpoint"],
        validate="1:1",
    )
    aligned.append(h.with_columns(pl.lit(method).alias("method")))
    for (ep,), g in h.group_by("endpoint"):
        for subset in ["all", "neighbors"]:
            x = g if subset == "all" else g.filter(pl.col("is_neighbor"))
            aux_rows.append(
                {
                    "method": method,
                    "head": ep,
                    "scope": "pred_vs_pic50_" + subset,
                    "n": x.height,
                    "RMSE": None,
                    "Spearman": float(spearmanr(x["pred_log2fc"], x["y_true"]).statistic),
                }
            )
pl.DataFrame(aux_rows).sort("scope", "head", "method").write_csv(out / "aux-metrics.csv")
pl.concat(aligned).write_parquet(out / "aligned.parquet")
write_json(
    out / "summary.json", {"config": cfg, "scope": "Untouched fold 0 only; no five-fold HPO claim"}
)
print(macro)
trials = json.loads((inp / cfg["trials"]).read_text())
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
axes[0].scatter([r["trial"] for r in trials], [r["value"] for r in trials], color="#356b99")
axes[0].axhline(trials[0]["value"], color="#b45d35", linestyle="--", label="Original")
axes[0].set(
    xlabel="Rustuna trial",
    ylabel="Best internal validation MSE (standardized)",
    title="16-trial search",
)
axes[0].legend()
for name, path in cfg["histories"].items():
    history = [json.loads(line) for line in (inp / path).read_text().splitlines()]
    for loss, style in [("train_mse_z", "--"), ("val_mse_z", "-")]:
        axes[1].plot(
            [r["epoch"] + 1 for r in history],
            [r[loss] for r in history],
            linestyle=style,
            color="#b45d35" if name == "original" else "#356b99",
            label=name + " " + ("train" if loss.startswith("train") else "validation"),
        )
axes[1].set(
    xlabel="Epoch", ylabel="Standardized MSE", title="Original versus selected learning curves"
)
axes[1].legend()
fig.savefig(out / "tuning.png", dpi=160)
plt.close(fig)
