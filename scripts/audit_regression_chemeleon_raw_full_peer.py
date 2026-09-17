"""Evaluate raw CheMeleon as a complete equal-weight ensemble member."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
values = ["y_true", "conf_low", "conf_high", "is_neighbor"]


def component(name: str, prediction_name: str, include_values: bool = False) -> pl.DataFrame:
    columns = [*keys, *values] if include_values else keys
    return (
        pl.read_parquet(inp / cfg[name])
        .sort(keys)
        .select(*columns, pl.col("prediction").alias(prediction_name))
    )


table = component("raw", "raw", include_values=True)
for name in ["pca", "frozen", "scores", "finetuned"]:
    table = table.join(component(name, name), on=keys, validate="1:1")
assert table.height == 6384 and table.select(keys).n_unique() == 6384

top2 = (
    pl.read_parquet(inp / cfg["cyp3a4_top2"])
    .filter(
        (pl.col("strategy") == "crossfold_topk")
        & (pl.col("weighting") == "weighted")
        & (pl.col("k") == 2)
    )
    .sort(keys)
    .select(*keys, pl.col("five_component_prediction").alias("current_3a4"))
)
assert top2.height == table.filter(pl.col("endpoint").str.starts_with("CYP3A4")).height
table = table.join(top2, on=keys, how="left", validate="m:1")

endpoint_components = {
    "CYP1A2_pIC50_direct_inhibition": ["frozen", "scores", "finetuned"],
    "CYP2C9_pIC50_direct_inhibition": ["scores", "finetuned", "pca"],
    "CYP2D6_pIC50_direct_inhibition": ["scores", "finetuned"],
}
parts = []
for endpoint, members in endpoint_components.items():
    dual_members = list(dict.fromkeys([*members, "raw", "pca"]))
    part = table.filter(pl.col("endpoint") == endpoint).with_columns(
        pl.mean_horizontal(*members).alias("current"),
        pl.mean_horizontal(*members, "raw").alias("raw_full_peer"),
        pl.mean_horizontal(*dual_members).alias("dual_view_full_peers"),
    )
    parts.append(part)
part = table.filter(pl.col("endpoint").str.starts_with("CYP3A4")).with_columns(
    pl.col("current_3a4").alias("current"),
    ((pl.col("current_3a4") * 5 + pl.col("raw")) / 6).alias("raw_full_peer"),
    ((pl.col("current_3a4") * 5 + pl.col("raw")) / 6).alias("dual_view_full_peers"),
)
parts.append(part)
aligned = pl.concat(parts).sort(keys)
methods = ["current", "raw_full_peer", "dual_view_full_peers"]
assert aligned.height == 6384 and aligned.select(keys).n_unique() == 6384
assert aligned.select(methods).select(pl.all().is_finite().all()).row(0) == (True,) * 3

spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
for (endpoint,), endpoint_table in aligned.group_by("endpoint", maintain_order=True):
    for fold in [-1, 0, 1, 2, 3, 4]:
        fold_table = endpoint_table if fold == -1 else endpoint_table.filter(pl.col("fold") == fold)
        y = fold_table["y_true"].to_numpy()
        for method in methods:
            prediction = fold_table[method].to_numpy()
            rows.append(
                {
                    "endpoint": endpoint,
                    "fold": fold,
                    "n": fold_table.height,
                    "method": method,
                    "ST_RAE": float(
                        scorer.rae_soft_threshold_absolute_error(
                            y,
                            prediction,
                            y_true_lower=fold_table["conf_low"].to_numpy(),
                            y_true_upper=fold_table["conf_high"].to_numpy(),
                        )
                    ),
                    "MAE": float(np.abs(y - prediction).mean()),
                    "Spearman": float(spearmanr(y, prediction).statistic),
                }
            )

metrics = pl.DataFrame(rows).sort("endpoint", "fold", "method")
pooled = metrics.filter(pl.col("fold") == -1)
macro = pooled.group_by("method").agg(pl.col("ST_RAE", "MAE", "Spearman").mean()).sort("ST_RAE")
fold_deltas = (
    metrics.filter(pl.col("fold") >= 0)
    .pivot(on="method", index=["endpoint", "fold"], values="ST_RAE")
    .with_columns(
        (pl.col("raw_full_peer") - pl.col("current")).alias("raw_peer_minus_current"),
        (pl.col("dual_view_full_peers") - pl.col("current")).alias("dual_peers_minus_current"),
    )
    .sort("endpoint", "fold")
)

aligned.write_parquet(out / "aligned-oof.parquet")
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
fold_deltas.write_csv(out / "fold-deltas.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": aligned.height,
            "methods": methods,
            "weight_search": False,
            "macro": macro.to_dicts(),
        },
        indent=2,
    )
)
print(macro)
print(pooled.sort("endpoint", "ST_RAE"))
print(fold_deltas)
