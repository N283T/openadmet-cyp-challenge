"""Compare paired parent-retention predictions on identical validation molecules."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.concat([pl.read_parquet(inp / p) for p in cfg["oof"]])
a = pl.read_parquet(inp / cfg["assignments"])
labels = pl.read_parquet(inp / cfg["labels"])
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
rows = []
boots = []
for ep in sorted(oof["endpoint"].unique()):
    ref = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .filter(~pl.col("is_anchor"))
        .sort("compound_id")
    )
    assert oof.filter(pl.col("endpoint") == ep).height == ref.height * 4
    for features in ["full", "mordred"]:
        maps = {}
        for mode in ["parent_excluded", "parent_retained"]:
            pred = oof.filter(
                (pl.col("endpoint") == ep)
                & (pl.col("features") == features)
                & (pl.col("mode") == mode)
            )
            f = ref.join(
                pred.select(
                    "compound_id", pl.col("fold").alias("pred_fold"), "y_true", "prediction"
                ),
                on="compound_id",
                validate="1:1",
                maintain_order="left",
            )
            assert (f["fold"] == f["pred_fold"]).all() and np.allclose(f["value"], f["y_true"])
            assert f["prediction"].is_finite().all()
            maps[mode] = f
            for subset in ["all_nonanchors", "neighbors"]:
                g = f if subset == "all_nonanchors" else f.filter(pl.col("is_neighbor"))
                for fold in [-1, 0, 1, 2, 3, 4]:
                    h = g if fold == -1 else g.filter(pl.col("fold") == fold)
                    y, p, lo, hi = [
                        h[c].to_numpy() for c in ["value", "prediction", "conf_low", "conf_high"]
                    ]
                    soft = np.maximum(lo - p, 0) + np.maximum(p - hi, 0)
                    den = np.maximum(lo - y.mean(), 0) + np.maximum(y.mean() - hi, 0)
                    rows.append(
                        dict(
                            endpoint=ep,
                            features=features,
                            mode=mode,
                            subset=subset,
                            fold=fold,
                            n=h.height,
                            ST_RAE=float(
                                official.rae_soft_threshold_absolute_error(
                                    y, p, y_true_lower=lo, y_true_upper=hi
                                )
                            ),
                            MAE=float(np.abs(y - p).mean()),
                            Spearman=float(spearmanr(y, p).statistic),
                            soft_MAE=float(soft.mean()),
                            denominator=float(den.mean()),
                            bias=float((p - y).mean()),
                            y_mean=float(y.mean()),
                            prediction_mean=float(p.mean()),
                        )
                    )
        f = maps["parent_retained"].filter(pl.col("is_neighbor"))
        ex = maps["parent_excluded"].filter(pl.col("is_neighbor"))
        assert f["compound_id"].to_list() == ex["compound_id"].to_list()
        y, lo, hi, p = [f[c].to_numpy() for c in ["value", "conf_low", "conf_high", "prediction"]]
        q = ex["prediction"].to_numpy()
        groups = f["group"].to_list()
        unique = sorted(set(groups))
        rng = np.random.default_rng(42)
        counts = rng.multinomial(len(unique), np.ones(len(unique)) / len(unique), size=1000)
        w = counts[:, [unique.index(g) for g in groups]]
        mu = w @ y / w.sum(axis=1)
        den = (
            w
            * (np.maximum(lo[None, :] - mu[:, None], 0) + np.maximum(mu[:, None] - hi[None, :], 0))
        ).sum(axis=1)
        delta = (
            w
            @ (
                np.maximum(lo - p, 0)
                + np.maximum(p - hi, 0)
                - np.maximum(lo - q, 0)
                - np.maximum(q - hi, 0)
            )
            / den
        )
        low, high = np.quantile(delta, [0.025, 0.975])
        boots.append(
            dict(endpoint=ep, features=features, delta_low=float(low), delta_high=float(high))
        )
m = pl.DataFrame(rows)
with psycopg.connect(**cfg["database"]) as conn:
    for r in rows:
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "parent_cv:" + r["subset"],
                    str(r["fold"]),
                    r["features"] + ":" + r["mode"] + ":" + metric,
                    r[metric],
                ),
            )
m.write_csv(out / "metrics.csv")
m.filter(pl.col("fold") == -1).group_by("features", "mode", "subset").agg(
    pl.col("ST_RAE", "MAE", "Spearman").mean()
).write_csv(out / "macro.csv")
pl.DataFrame(boots).write_csv(out / "paired-bootstrap.csv")
oof.write_parquet(out / "oof.parquet")
write_json(
    out / "summary.json",
    dict(
        config=cfg,
        key_label_fold_checks=True,
        comparison=(
            "Matched molecules, nested auxiliaries and training-only PCA/imputation. "
            "Fixed LightGBM; bootstrap conditional on fitted models."
        ),
    ),
)
print(m.filter((pl.col("fold") == -1) & (pl.col("subset") == "neighbors")))
