"""Evaluate direct and TabPFN fine-tuned variants with matched fixed blends."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def p(n):
    return next(inp.glob("*-" + n))


keys = ["compound_id", "endpoint", "fold"]
base = pl.read_parquet(p("aligned-oof.parquet")).sort(keys)
new = pl.read_parquet(p("finetune-oof.parquet"))
new_methods = ["admetica_ft_direct", "admetica_ft_emb_tabpfn"]
for method in new_methods:
    h = new.filter(pl.col("method") == method).sort(keys)
    assert h.height == 6384 and not h.select(keys).is_duplicated().any()
    for c in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert h[c].equals(base[c]), (method, c)
    base = base.with_columns(h["prediction"].alias(method))
methods = [
    "admetica_emb_tabpfn",
    "mean_fp_log2fc_admetica",
    "mean__fp_mordred_lgbm__pred_log2fc_tabpfn",
    "current_full_tabpfn",
] + new_methods
for method in new_methods:
    for pool, names in [
        ("triple", ["fp_mordred_lgbm", "pred_log2fc_tabpfn", method]),
        ("four", ["fp_mordred_lgbm", "pred_log2fc_tabpfn", "admetica_emb_tabpfn", method]),
    ]:
        name = pool + "_" + method
        base = base.with_columns(pl.mean_horizontal(*names).alias(name))
        methods.append(name)
base.write_parquet(out / "aligned-oof.parquet")
spec = importlib.util.spec_from_file_location("official", p("custom_scoring_functions.py"))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows, correlations = [], []
for (ep,), group in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        selected = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = selected if fold == -1 else selected.filter(pl.col("fold") == fold)
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
            if fold == -1:
                for method in new_methods:
                    for other in ["fp_mordred_lgbm", "pred_log2fc_tabpfn", "admetica_emb_tabpfn"]:
                        x, z = h[method].to_numpy(), h[other].to_numpy()
                        correlations.append(
                            dict(
                                endpoint=ep,
                                subset=subset,
                                method=method,
                                other=other,
                                pearson=float(np.corrcoef(x, z)[0, 1]),
                                residual_pearson=float(np.corrcoef(x - y, z - y)[0, 1]),
                            )
                        )
m = pl.DataFrame(rows)
m.write_csv(out / "metrics.csv")
pl.DataFrame(correlations).write_csv(out / "correlations.csv")
macro = (
    m.filter(pl.col("fold") == -1)
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
                    "admetica_ft:" + r["subset"],
                    str(r["fold"]),
                    r["method"] + ":" + name,
                    r[name],
                ),
            )
(out / "summary.json").write_text(
    json.dumps(dict(config=cfg, matched_rows_per_method=6384, weight_search=False), indent=2)
)
print(macro)
