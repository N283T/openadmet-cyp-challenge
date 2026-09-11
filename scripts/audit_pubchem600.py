"""Matched PubChem600 embedding OOF scores, correlations, and fixed equal blends."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def p(name):
    return next(inp.glob("*-" + name))


base = pl.read_parquet(p("aligned-oof.parquet"))
new = pl.read_parquet(p("pubchem600-oof.parquet"))
keys = ["compound_id", "endpoint", "fold"]
base, new = base.sort(keys), new.sort(keys)
assert new.height == 6384 and not new.select(keys).is_duplicated().any()
for c in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert base[c].equals(new[c]), c
base = base.with_columns(new["prediction"].alias("pubchem600_tabpfn"))
methods = [
    "fp_mordred_lgbm",
    "emb_tabpfn",
    "pred_log2fc_tabpfn",
    "pubchem600_tabpfn",
    "current_full_tabpfn",
    "current_full_tabicl",
    "mean__fp_mordred_lgbm__pred_log2fc_tabpfn",
    "admetica_emb_tabpfn",
    "mean_fp_log2fc_admetica",
]
for old in methods[:3]:
    name = "mean_pubchem600_" + old
    base = base.with_columns(((pl.col(old) + pl.col("pubchem600_tabpfn")) / 2).alias(name))
    methods.append(name)
base = base.with_columns(
    pl.mean_horizontal("fp_mordred_lgbm", "pred_log2fc_tabpfn", "pubchem600_tabpfn").alias(
        "mean_fp_log2fc_pubchem600"
    )
)
methods.append("mean_fp_log2fc_pubchem600")
base = base.with_columns(
    pl.mean_horizontal(
        "fp_mordred_lgbm", "pred_log2fc_tabpfn", "admetica_emb_tabpfn", "pubchem600_tabpfn"
    ).alias("mean_four_pubchem600")
)
methods.append("mean_four_pubchem600")
base.write_parquet(out / "aligned-oof.parquet")
spec = importlib.util.spec_from_file_location("official", p("custom_scoring_functions.py"))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows, corr = [], []
for (ep,), g in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        group = g if subset == "all" else g.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = group if fold == -1 else group.filter(pl.col("fold") == fold)
            y = h["y_true"].to_numpy()
            for method in methods:
                pred = h[method].to_numpy()
                rows.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        n=h.height,
                        method=method,
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                y,
                                pred,
                                y_true_lower=h["conf_low"].to_numpy(),
                                y_true_upper=h["conf_high"].to_numpy(),
                            )
                        ),
                        MAE=float(np.abs(y - pred).mean()),
                        Spearman=float(spearmanr(y, pred).statistic),
                    )
                )
            for old in methods[:3]:
                x, z = h[old].to_numpy(), h["pubchem600_tabpfn"].to_numpy()
                corr.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        other=old,
                        pearson=float(np.corrcoef(x, z)[0, 1]),
                        spearman=float(spearmanr(x, z).statistic),
                        residual_pearson=float(np.corrcoef(x - y, z - y)[0, 1]),
                    )
                )
metrics = pl.DataFrame(rows)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(corr).write_csv(out / "correlations.csv")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
macro.write_csv(out / "macro.csv")
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
with psycopg.connect(**cfg["database"]) as conn:
    for r in rows:
        for name in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "pubchem600:" + r["subset"],
                    str(r["fold"]),
                    r["method"] + ":" + name,
                    r[name],
                ),
            )
(out / "summary.json").write_text(
    json.dumps(
        dict(config=cfg, matched_oof=6384, models=20, ensemble_weight_search=False), indent=2
    )
)
print(macro)
