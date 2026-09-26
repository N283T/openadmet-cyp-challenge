"""Compare a recorded Live/Interim CYP leaderboard snapshot without hidden labels."""

import hashlib
import json
import os
from pathlib import Path

import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.tracks import ENDPOINTS


def main():
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    subject = config["subject"]
    frames = {}
    for entry in config["entries"]:
        key = (entry["board"], entry["track"], entry["endpoint"])
        path = input_dir / entry["alias"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Input hash changed: {entry['alias']}")
        frame = pl.read_csv(path)
        if frame.is_empty() or frame["Username"].n_unique() != frame.height:
            raise ValueError(f"Empty or duplicate-user leaderboard: {key}")
        frames[key] = frame

    expected = {
        (board, track, endpoint)
        for board in ("live", "interim")
        for track in ENDPOINTS
        for endpoint in ("MA", *ENDPOINTS[track])
    }
    if set(frames) != expected:
        raise ValueError(f"Missing or unexpected boards: {expected ^ set(frames)}")

    own = []
    stability = []
    head_to_head = []
    for track, endpoints in ENDPOINTS.items():
        for endpoint in ("MA", *endpoints):
            metric = ("MA-" if endpoint == "MA" else "") + (
                "ST-RAE" if track == "regression" else "MCC"
            )
            live = frames["live", track, endpoint]
            interim = frames["interim", track, endpoint]
            for board, frame in (("live", live), ("interim", interim)):
                subject_row = frame.filter(pl.col("Username") == subject)
                if subject_row.height != 1:
                    raise ValueError(f"Subject missing or duplicated: {board}, {track}, {endpoint}")
                row = subject_row.row(0, named=True)
                leader = frame.row(0, named=True)
                own.append(
                    {
                        "board": board,
                        "track": track,
                        "endpoint": endpoint,
                        "metric": metric,
                        "participants": frame.height,
                        "rank": row["Rank"],
                        "score": row[metric],
                        "leader": leader["Username"],
                        "leader_score": leader[metric],
                        "gap_to_leader": (
                            row[metric] - leader[metric]
                            if track == "regression"
                            else leader[metric] - row[metric]
                        ),
                        "tier": row.get("Significance (tiers)"),
                    }
                )
            paired = live.select("Username", pl.col("Rank").alias("live_rank")).join(
                interim.select("Username", pl.col("Rank").alias("interim_rank")),
                on="Username",
                how="inner",
            )
            correlation = spearmanr(paired["live_rank"], paired["interim_rank"]).statistic
            top_live = set(live.head(10)["Username"])
            top_interim = set(interim.head(10)["Username"])
            stability.append(
                {
                    "track": track,
                    "endpoint": endpoint,
                    "common_participants": paired.height,
                    "rank_spearman": correlation,
                    "top10_overlap": len(top_live & top_interim),
                    "live_count": live.height,
                    "interim_count": interim.height,
                }
            )
        comparator = frames["interim", track, "MA"].row(0, named=True)["Username"]
        for endpoint in ("MA", *endpoints):
            frame = frames["interim", track, endpoint]
            metric = ("MA-" if endpoint == "MA" else "") + (
                "ST-RAE" if track == "regression" else "MCC"
            )
            scores = frame.filter(pl.col("Username").is_in([subject, comparator]))
            if scores.height != (1 if subject == comparator else 2):
                raise ValueError(f"Comparator missing: {track}, {endpoint}")
            own_score = scores.filter(pl.col("Username") == subject)[metric].item()
            comparator_score = scores.filter(pl.col("Username") == comparator)[metric].item()
            head_to_head.append(
                {
                    "track": track,
                    "endpoint": endpoint,
                    "subject": subject,
                    "comparator": comparator,
                    "subject_score": own_score,
                    "comparator_score": comparator_score,
                    "subject_minus_comparator": own_score - comparator_score,
                }
            )

    pl.DataFrame(own).sort("track", "endpoint", "board").write_csv(output_dir / "own.csv")
    pl.DataFrame(stability).sort("track", "endpoint").write_csv(
        output_dir / "rank_stability.csv"
    )
    pl.DataFrame(head_to_head).sort("track", "endpoint").write_csv(
        output_dir / "head_to_head.csv"
    )
    (output_dir / "summary.json").write_text(
        json.dumps(
            {
                "subject": subject,
                "source_manifest": config["source_manifest"],
                "boards": len(frames),
                "scope": "Descriptive leaderboard comparison; no hidden labels or predictions",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
