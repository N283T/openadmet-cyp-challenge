"""Describe each CYP direct-inhibition endpoint on frozen Live/Interim boards."""

import hashlib
import json
import os
from pathlib import Path

import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.tracks import ENDPOINTS


def metric_rank(frame, value, column, higher_is_better=False):
    better = pl.col(column) > value if higher_is_better else pl.col(column) < value
    return frame.filter(better).height + 1


def main():
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    frames = {}
    for entry in config["entries"]:
        if entry["track"] != "regression":
            raise ValueError("Regression inputs only")
        path = input_dir / entry["alias"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Input hash changed: {entry['alias']}")
        frame = pl.read_csv(path)
        if frame.is_empty() or frame["Username"].n_unique() != frame.height:
            raise ValueError(f"Empty or duplicate-user leaderboard: {entry['alias']}")
        frames[entry["board"], entry["endpoint"]] = frame

    endpoints = ENDPOINTS["regression"]
    expected = {(board, ep) for board in ("live", "interim") for ep in ("MA", *endpoints)}
    if set(frames) != expected:
        raise ValueError(f"Coverage mismatch: {expected ^ set(frames)}")

    subject = config["subject"]
    target_rows = []
    top_rows = []
    for endpoint in endpoints:
        live = frames["live", endpoint]
        interim = frames["interim", endpoint]
        mine = interim.filter(pl.col("Username") == subject).row(0, named=True)
        mine_live = live.filter(pl.col("Username") == subject).row(0, named=True)
        leader = interim.row(0, named=True)
        paired = live.select("Username", pl.col("Rank").alias("live_rank")).join(
            interim.select("Username", pl.col("Rank").alias("interim_rank")),
            on="Username",
            how="inner",
        )
        target_rows.append(
            {
                "endpoint": endpoint,
                "live_rank": mine_live["Rank"],
                "interim_rank": mine["Rank"],
                "interim_entrants": interim.height,
                "interim_rank_percentile": mine["Rank"] / interim.height,
                "live_st_rae": mine_live["ST-RAE"],
                "interim_st_rae": mine["ST-RAE"],
                "leader": leader["Username"],
                "leader_st_rae": leader["ST-RAE"],
                "gap_to_leader": mine["ST-RAE"] - leader["ST-RAE"],
                "interim_median_st_rae": interim["ST-RAE"].median(),
                "interim_q10_st_rae": interim["ST-RAE"].quantile(0.1),
                "interim_mae": mine["MAE"],
                "mae_rank": metric_rank(interim, mine["MAE"], "MAE"),
                "interim_spearman": mine["Spearman's ρ"],
                "spearman_rank": metric_rank(
                    interim, mine["Spearman's ρ"], "Spearman's ρ", True
                ),
                "leader_mae": leader["MAE"],
                "leader_spearman": leader["Spearman's ρ"],
                "common_entrants": paired.height,
                "live_interim_rank_spearman": spearmanr(
                    paired["live_rank"], paired["interim_rank"]
                ).statistic,
            }
        )
        for row in interim.head(10).to_dicts():
            top_rows.append(
                {
                    "endpoint": endpoint,
                    "rank": row["Rank"],
                    "username": row["Username"],
                    "st_rae": row["ST-RAE"],
                    "mae": row["MAE"],
                    "spearman": row["Spearman's ρ"],
                }
            )

    pl.DataFrame(target_rows).write_csv(output_dir / "targets.csv")
    pl.DataFrame(top_rows).write_csv(output_dir / "top10.csv")
    (output_dir / "summary.json").write_text(
        json.dumps(
            {
                "subject": subject,
                "source_manifest": config["source_manifest"],
                "interpretation": "Descriptive board statistics; no hidden labels or paired errors",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
