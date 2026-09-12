"""Audit matched external readouts and prespecified equal ensemble additions."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(name):
    matches = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
    assert len(matches) == 1, (name, matches)
    return matches[0]


keys = ["compound_id", "endpoint", "fold"]
base = pl.read_parquet(path("aligned-oof.parquet")).sort(keys)
new = pl.read_parquet(path(cfg.get("prediction_file", "external-oof.parquet"))).sort(keys)
assert new.height == base.height == 6384
assert not new.select(keys).is_duplicated().any()
for c in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert base[c].equals(new[c]), c
assert new["method"].n_unique() == 1
method = new["method"][0]
base = base.with_columns(new["prediction"].alias(method))
assert np.isfinite(new["prediction"].to_numpy()).all()
controls = [
    "fp_mordred_lgbm",
    "emb_tabpfn",
    "pred_log2fc_tabpfn",
    "admetica_emb_tabpfn",
    "pubchem600_tabpfn",
    "mean_fp_log2fc_admetica",
    "mean_four_pubchem600",
]
methods = controls + [method]
for suffix, old in [
    ("triple", ["fp_mordred_lgbm", "pred_log2fc_tabpfn"]),
    ("four", ["fp_mordred_lgbm", "pred_log2fc_tabpfn", "admetica_emb_tabpfn"]),
    ("five", ["fp_mordred_lgbm", "pred_log2fc_tabpfn", "admetica_emb_tabpfn", "pubchem600_tabpfn"]),
]:
    name = method + "__" + suffix
    base = base.with_columns(pl.mean_horizontal(*old, method).alias(name))
    methods.append(name)
spec = importlib.util.spec_from_file_location("official", path("custom_scoring_functions.py"))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows, correlations = [], []
for (ep,), g in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        group = g if subset == "all" else g.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = group if fold == -1 else group.filter(pl.col("fold") == fold)
            y = h["y_true"].to_numpy()
            for name in methods:
                p = h[name].to_numpy()
                rows.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        n=h.height,
                        method=name,
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                y,
                                p,
                                y_true_lower=h["conf_low"].to_numpy(),
                                y_true_upper=h["conf_high"].to_numpy(),
                            )
                        ),
                        MAE=float(np.abs(y - p).mean()),
                        Spearman=float(spearmanr(y, p).statistic),
                    )
                )
            for old in controls[:5]:
                x, z = h[old].to_numpy(), h[method].to_numpy()
                correlations.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        other=old,
                        pearson=float(np.corrcoef(x, z)[0, 1]),
                        residual_pearson=float(np.corrcoef(x - y, z - y)[0, 1]),
                    )
                )
metrics = pl.DataFrame(rows)
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
pl.DataFrame(correlations).write_csv(out / "correlations.csv")
base.write_parquet(out / "aligned-oof.parquet")
with psycopg.connect(**cfg["database"]) as conn:
    for r in rows:
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "external:" + r["subset"],
                    str(r["fold"]),
                    r["method"] + ":" + metric,
                    r[metric],
                ),
            )
(out / "summary.json").write_text(
    json.dumps(dict(config=cfg, matched_oof=new.height, weight_search=False), indent=2)
)
print(macro)
