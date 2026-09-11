"""Measure auxiliary OOF prediction association with direct pIC50."""

import os
from pathlib import Path

import polars as pl
from scipy.stats import pearsonr, spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
frame = pl.read_parquet(inp / "000-aligned-residuals.parquet")
assert frame.select("compound_id", "endpoint", "fold").is_duplicated().sum() == 0
rows = []
for (endpoint,), group in frame.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        selected = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        rows.append(
            {
                "endpoint": endpoint,
                "subset": subset,
                "n": selected.height,
                "pred_log2fc_pic50_spearman": float(
                    spearmanr(selected["predicted_primary"], selected["y_true"]).statistic
                ),
                "pred_log2fc_pic50_pearson": float(
                    pearsonr(selected["predicted_primary"], selected["y_true"]).statistic
                ),
                "observed_log2fc_pic50_spearman": float(
                    spearmanr(selected["observed_primary"], selected["y_true"]).statistic
                ),
            }
        )
result = pl.DataFrame(rows).sort("endpoint", "subset")
result.write_csv(out / "correlations.csv")
print(result)
