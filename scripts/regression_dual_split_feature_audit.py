"""Compare matched feature candidates across primary and secondary regression splits."""

import json
import os
from pathlib import Path

import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

primary = pl.read_csv(inp / "primary-metrics.csv").select("method", "endpoint", "ST_RAE")
secondary = pl.read_csv(inp / "secondary-metrics.csv").select("method", "endpoint", "ST_RAE")
common = sorted(set(primary["method"]) & set(secondary["method"]))
primary = primary.filter(pl.col("method").is_in(common))
secondary = secondary.filter(pl.col("method").is_in(common))

joined = primary.rename({"ST_RAE": "primary_ST_RAE"}).join(
    secondary.rename({"ST_RAE": "secondary_ST_RAE"}),
    on=["method", "endpoint"],
    validate="1:1",
)
if joined.height != len(common) * primary["endpoint"].n_unique():
    raise ValueError("Matched candidate grid is incomplete")

comparison_rows = []
summary_rows = []
for endpoint in sorted(joined["endpoint"].unique()):
    frame = joined.filter(pl.col("endpoint") == endpoint)
    control = frame.filter(pl.col("method") == cfg["control"])
    if control.height != 1:
        raise ValueError(f"Missing control for {endpoint}")
    primary_control = control["primary_ST_RAE"][0]
    secondary_control = control["secondary_ST_RAE"][0]
    frame = frame.with_columns(
        pl.col("primary_ST_RAE").rank("ordinal").cast(pl.Int64).alias("primary_rank"),
        pl.col("secondary_ST_RAE").rank("ordinal").cast(pl.Int64).alias("secondary_rank"),
        (pl.col("primary_ST_RAE") - primary_control).alias("primary_delta_vs_control"),
        (pl.col("secondary_ST_RAE") - secondary_control).alias("secondary_delta_vs_control"),
    ).with_columns(
        (
            (pl.col("primary_delta_vs_control") < 0) & (pl.col("secondary_delta_vs_control") < 0)
        ).alias("improves_control_both"),
        (pl.col("primary_rank") - pl.col("secondary_rank")).abs().alias("absolute_rank_shift"),
    )
    comparison_rows.append(frame)
    p_winner = frame.sort("primary_ST_RAE").row(0, named=True)
    s_winner = frame.sort("secondary_ST_RAE").row(0, named=True)
    summary_rows.append(
        {
            "endpoint": endpoint,
            "candidates": frame.height,
            "score_rank_spearman": float(
                spearmanr(frame["primary_ST_RAE"], frame["secondary_ST_RAE"]).statistic
            ),
            "primary_winner": p_winner["method"],
            "primary_winner_primary_ST_RAE": p_winner["primary_ST_RAE"],
            "primary_winner_secondary_ST_RAE": p_winner["secondary_ST_RAE"],
            "primary_winner_secondary_rank": p_winner["secondary_rank"],
            "secondary_winner": s_winner["method"],
            "secondary_winner_primary_ST_RAE": s_winner["primary_ST_RAE"],
            "secondary_winner_secondary_ST_RAE": s_winner["secondary_ST_RAE"],
            "secondary_winner_primary_rank": s_winner["primary_rank"],
            "improves_control_both": int(frame["improves_control_both"].sum()),
        }
    )

comparison = pl.concat(comparison_rows).sort("endpoint", "primary_rank")
summary = pl.DataFrame(summary_rows).sort("endpoint")
comparison.write_csv(out / "candidate-comparison.csv")
summary.write_csv(out / "endpoint-summary.csv")
write_json(
    out / "summary.json",
    {
        "primary_split": cfg["primary_split"],
        "secondary_split": cfg["secondary_split"],
        "control": cfg["control"],
        "policy": cfg["policy"],
        "common_candidates": len(common),
        "excluded_primary_only": sorted(
            set(pl.read_csv(inp / "primary-metrics.csv")["method"]) - set(common)
        ),
        "endpoints": summary.to_dicts(),
    },
)
print(summary, flush=True)
