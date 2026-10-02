"""Is the CYP3A4 TDI Live rank 1 distinguishable from the top tier?

Compares the frozen September 26 Live and Interim CYP3A4_is_TDI leaderboards on the
entrants common to both boards, expresses N283T's Live lead and Interim deficit in
units of the reported bootstrap standard deviations, infers per-entrant scored
positive counts and call fractions from the published aggregate metrics, and lists
the N283T A4 Live history across the saved snapshots.

No model is fitted and no prediction file is read; this run consumes only public
leaderboard CSVs.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

METRICS = [
    "MCC",
    "MCC std",
    "Accuracy",
    "Accuracy std",
    "Precision",
    "Precision std",
    "Recall",
    "Recall std",
    "F1 Score",
    "F1 Score std",
]


def load_board(path):
    frame = pl.read_csv(path, infer_schema_length=0)
    prefix = "MA-" if "MA-MCC" in frame.columns else ""
    renames = {f"{prefix}{m}": m for m in METRICS if f"{prefix}{m}" in frame.columns}
    frame = frame.rename(renames)
    return frame.select(
        pl.col("Rank").cast(pl.Int64).alias("rank"),
        pl.col("Username").alias("username"),
        pl.col("Submitted").alias("submitted"),
        *[pl.col(m).cast(pl.Float64).alias(m.lower().replace(" ", "_")) for m in METRICS],
    )


def infer_positives(recall, recall_std):
    """Bootstrap std of a proportion over P scored positives implies P = r(1-r)/std^2."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(recall_std > 0, recall * (1.0 - recall) / recall_std**2, np.nan)


def infer_operating_point(accuracy, precision, recall):
    """Scored prevalence and positive-call fraction from the published aggregates."""
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 1.0 - 2.0 * recall + recall / precision
        prevalence = np.where(np.abs(denom) > 1e-9, (1.0 - accuracy) / denom, np.nan)
        call_fraction = prevalence * recall / precision
    return prevalence, call_fraction


def operating_table(board, top_n, board_name):
    frame = board.sort("rank").head(top_n)
    recall = frame["recall"].to_numpy()
    positives = infer_positives(recall, frame["recall_std"].to_numpy())
    prevalence, call_fraction = infer_operating_point(
        frame["accuracy"].to_numpy(), frame["precision"].to_numpy(), recall
    )
    scored = np.where(prevalence > 0, positives / prevalence, np.nan)
    return frame.with_columns(
        pl.lit(board_name).alias("board"),
        pl.Series("inferred_scored_positives", positives),
        pl.Series("inferred_prevalence", prevalence),
        pl.Series("inferred_call_fraction", call_fraction),
        pl.Series("implied_scored_molecules", scored),
    ).select(
        "board",
        "rank",
        "username",
        "mcc",
        "mcc_std",
        "accuracy",
        "precision",
        "recall",
        "f1_score",
        "inferred_scored_positives",
        "inferred_prevalence",
        "inferred_call_fraction",
        "implied_scored_molecules",
    )


def spearman(a, b):
    return float(pl.DataFrame({"a": a, "b": b}).select(pl.corr("a", "b", method="spearman")).item())


def pearson(a, b):
    return float(pl.DataFrame({"a": a, "b": b}).select(pl.corr("a", "b")).item())


def board_comparison(live, interim, name, username, top_k=10):
    common = (
        live.select("username", pl.col("mcc").alias("live_mcc"), pl.col("rank").alias("live_rank"))
        .join(
            interim.select(
                "username",
                pl.col("mcc").alias("interim_mcc"),
                pl.col("rank").alias("interim_rank"),
            ),
            on="username",
            how="inner",
        )
        .sort("live_rank")
    )
    live_top = set(live.sort("rank").head(top_k)["username"].to_list())
    interim_top = set(interim.sort("rank").head(top_k)["username"].to_list())
    result = dict(
        endpoint=name,
        live_entrants=live.height,
        interim_entrants=interim.height,
        common_entrants=common.height,
        mcc_spearman=spearman(common["live_mcc"], common["interim_mcc"]),
        mcc_pearson=pearson(common["live_mcc"], common["interim_mcc"]),
        rank_spearman=spearman(common["live_rank"], common["interim_rank"]),
        top10_overlap=len(live_top & interim_top),
        mean_live_minus_interim=float(
            (common["live_mcc"] - common["interim_mcc"]).mean()  # type: ignore[operator]
        ),
        sd_live_minus_interim=float(
            (common["live_mcc"] - common["interim_mcc"]).std()  # type: ignore[operator]
        ),
    )
    row = common.filter(pl.col("username") == username)
    if row.height == 1:
        result.update(
            user_live_mcc=float(row["live_mcc"][0]),
            user_live_rank=int(row["live_rank"][0]),
            user_interim_mcc=float(row["interim_mcc"][0]),
            user_interim_rank=int(row["interim_rank"][0]),
        )
    return result, common


def separation(board, username, board_name):
    """N283T versus the board leader and its nearest rivals, in bootstrap std units."""
    ordered = board.sort("rank")
    mine = ordered.filter(pl.col("username") == username)
    assert mine.height == 1
    my_mcc, my_std = float(mine["mcc"][0]), float(mine["mcc_std"][0])
    my_rank = int(mine["rank"][0])
    leader = ordered.head(1)
    rows = []
    for label, other in [
        ("leader", leader),
        ("neighbour", ordered.filter(pl.col("rank") == (2 if my_rank == 1 else my_rank - 1))),
    ]:
        if other.height != 1:
            continue
        other_mcc, other_std = float(other["mcc"][0]), float(other["mcc_std"][0])
        diff = my_mcc - other_mcc
        pooled = float(np.hypot(my_std, other_std))
        rows.append(
            dict(
                board=board_name,
                relation=label,
                username=other["username"][0],
                rank=int(other["rank"][0]),
                mcc=other_mcc,
                mcc_std=other_std,
                user_mcc=my_mcc,
                user_mcc_std=my_std,
                user_rank=my_rank,
                delta=diff,
                delta_in_user_std=diff / my_std if my_std else None,
                delta_in_other_std=diff / other_std if other_std else None,
                delta_in_pooled_std=diff / pooled if pooled else None,
            )
        )
    within = {}
    for k in (1.0, 2.0):
        pooled = np.hypot(my_std, ordered["mcc_std"].to_numpy())
        close = np.abs(ordered["mcc"].to_numpy() - my_mcc) <= k * pooled
        within[f"entrants_within_{k:g}_pooled_std"] = int(close.sum())
        above = (ordered["mcc"].to_numpy() > my_mcc - k * pooled) & (
            ordered["username"].to_numpy() != username
        )
        within[f"entrants_above_minus_{k:g}_pooled_std"] = int(above.sum())
    return rows, within


def main():
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs, user = cfg["inputs"], cfg["username"]
    boards = {
        key: load_board(source / inputs[key])
        for key in ["live_a4", "live_d6", "live_ma", "interim_a4", "interim_d6", "interim_ma"]
    }
    manifest = json.loads((source / inputs["snapshot_manifest"]).read_text())

    # Board identity checks: both A4 boards must carry the expected N283T values.
    identity = {}
    for key, board in boards.items():
        row = board.filter(pl.col("username") == user)
        identity[key] = dict(
            entrants=board.height,
            user_rank=int(row["rank"][0]) if row.height else None,
            user_mcc=float(row["mcc"][0]) if row.height else None,
            user_submitted=row["submitted"][0] if row.height else None,
        )
    assert abs(identity["live_a4"]["user_mcc"] - 0.529911) < 1e-5
    assert abs(identity["interim_a4"]["user_mcc"] - 0.502251) < 1e-5
    print(
        json.dumps(dict(identity=identity, manifest_captured=manifest["captured_on"])), flush=True
    )

    comparisons, scatters = [], []
    for name, live_key, interim_key in [
        ("CYP3A4_is_TDI", "live_a4", "interim_a4"),
        ("CYP2D6_is_TDI", "live_d6", "interim_d6"),
        ("MA", "live_ma", "interim_ma"),
    ]:
        result, common = board_comparison(boards[live_key], boards[interim_key], name, user)
        comparisons.append(result)
        scatters.append(common.with_columns(pl.lit(name).alias("endpoint")))
    pl.concat(scatters).write_csv(out / "live_interim_scatter.csv")

    separations, neighbourhood = [], {}
    for board_name, key in [("live_a4", "live_a4"), ("interim_a4", "interim_a4")]:
        rows, within = separation(boards[key], user, board_name)
        separations.extend(rows)
        neighbourhood[board_name] = within
    pl.DataFrame(separations).write_csv(out / "separation.csv")

    operating = pl.concat(
        [
            operating_table(boards["live_a4"], cfg["top_n_operating_points"], "live_a4"),
            operating_table(boards["interim_a4"], cfg["top_n_operating_points"], "interim_a4"),
        ]
    )
    operating.write_csv(out / "operating_points.csv")

    # N283T A4 Live history across the saved snapshots.
    history_rows = []
    for submitted, alias in inputs["history"]:
        board = load_board(source / alias)
        row = board.filter(pl.col("username") == user)
        assert row.height == 1, alias
        assert row["submitted"][0] == submitted, (alias, row["submitted"][0], submitted)
        mcc = float(row["mcc"][0])
        assert any(abs(mcc - k) < 5e-4 for k in cfg["known_live_a4_mcc"]), (alias, mcc)
        history_rows.append(
            dict(
                submitted=submitted,
                alias=alias,
                entrants=board.height,
                rank=int(row["rank"][0]),
                mcc=mcc,
                mcc_std=float(row["mcc_std"][0]),
                accuracy=float(row["accuracy"][0]),
                precision=float(row["precision"][0]),
                recall=float(row["recall"][0]),
                f1=float(row["f1_score"][0]),
                recall_std=float(row["recall_std"][0]),
                leader_mcc=float(board.sort("rank")["mcc"][0]),
            )
        )
    history = pl.DataFrame(history_rows).with_columns(
        pl.col("mcc").diff().alias("mcc_delta"),
        (pl.col("mcc").diff() / pl.col("mcc_std")).alias("mcc_delta_in_std"),
    )
    recall = history["recall"].to_numpy()
    positives, call_fraction = infer_operating_point(
        history["accuracy"].to_numpy(), history["precision"].to_numpy(), recall
    )
    history = history.with_columns(
        pl.Series("inferred_prevalence", positives),
        pl.Series("inferred_call_fraction", call_fraction),
        pl.Series(
            "inferred_scored_positives", infer_positives(recall, history["recall_std"].to_numpy())
        ),
    )
    history.write_csv(out / "a4_live_history.csv")

    summary = dict(
        config=cfg,
        identity=identity,
        snapshot_captured_on=manifest["captured_on"],
        comparisons=comparisons,
        separation=separations,
        neighbourhood=neighbourhood,
        operating_points=operating.to_dicts(),
        history=history.to_dicts(),
        caveats=[
            "Published metrics are separately averaged bootstrap estimates, so inferred "
            "positive counts, prevalences and call fractions are approximations",
            "The reported std is a bootstrap resampling spread for one fixed prediction "
            "vector, not a between-entrant sampling error",
            "Live scores one half of the blinded test set; Interim scores the full set once",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(dict(comparisons=comparisons, neighbourhood=neighbourhood)), flush=True)


if __name__ == "__main__":
    main()
