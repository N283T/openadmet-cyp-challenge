"""Pre-register the next D6 probe and its prediction before any upload exists.

Parts 1 to 3 of the pseudo-negative validation leave one open, uploadable question: the
positive-call count. The rate grid, once extended past the 10% boundary that truncated
every earlier selection, has an interior optimum well above the deployed 49 calls, and the
two split families disagree about where it is. This script fixes the candidate call vectors,
counts exactly how many calls each one changes against the deployed October 3 payload, and
writes down the public-metric consequence of each as an explicit function of the marginal
precision of the added ranks.

The reconstructed Live state of the deployed payload is taken as given from run
20261003T012000305555Z-bca66716: 374 scored rows of 750, 26 scored positives, and for the
October 3 vector 28 scored calls with 11 true positives. Those are reconstructions from the
published bootstrap aggregates, not disclosed labels.

Nothing is fitted. No candidate, no prepared submission, no upload, and no API contact. The
output is a prediction to be checked later, which is what makes the next transition an
out-of-sample test rather than another retrospective fit.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl

ENDPOINT = "CYP2D6_is_TDI"


def mcc_from_counts(tp, fp, fn, tn):
    tp, fp, fn, tn = (float(v) for v in (tp, fp, fn, tn))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denom) if denom > 0 else 0.0


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def live_mcc(state, added_scored, marginal_precision):
    """Public MCC after adding `added_scored` scored calls at a given marginal precision."""
    gained = added_scored * marginal_precision
    tp = state["true_positives"] + gained
    scored_calls = state["scored_calls"] + added_scored
    fp = scored_calls - tp
    fn = state["scored_positives"] - tp
    tn = state["scored_rows"] - scored_calls - fn
    return mcc_from_counts(tp, fp, fn, tn)


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = cfg["inputs"]
    state = cfg["live_state"]
    deployed = pl.read_parquet(source / inputs["payload_oct03"]).select(
        pl.col("Molecule_Name").alias("molecule_name"), pl.col(ENDPOINT).alias("call")
    )
    assert deployed.height == 750
    deployed_calls = set(deployed.filter(pl.col("call"))["molecule_name"].to_list())
    assert len(deployed_calls) == cfg["deployed_calls"]

    legacy = pl.read_parquet(source / inputs["sep15_blind_scores"])
    legacy = legacy.filter(pl.col("model") == cfg["legacy_model"]).select("compound_id", "score")
    assert legacy.height == 750
    pool = pl.read_parquet(source / inputs["pool_blind_scores"])
    names = (
        pl.read_parquet(source / inputs["test_mapping"])
        .select("compound_id", "molecule_name")
        .unique()
    )
    pool = (
        pool.filter(pl.col("model") == "pool_lgbm")
        .select("compound_id", pl.col("score").alias("pool_score"))
        .join(legacy.rename({"score": "legacy_score"}), on="compound_id", validate="1:1")
        .join(names, on="compound_id", how="left", validate="1:1")
    )
    assert pool.height == 750 and pool["molecule_name"].null_count() == 0
    l4_score = (
        percentile_rank(pool["legacy_score"].to_numpy())
        + percentile_rank(pool["pool_score"].to_numpy())
    ) / 2
    combined = pool.with_columns(pl.Series("l4_score", l4_score))

    # The reconstructed L4 ordering must reproduce the deployed 49 calls exactly.
    ordered = combined.sort(["l4_score", "compound_id"], descending=[True, False])
    top_deployed = set(ordered["molecule_name"].to_list()[: cfg["deployed_calls"]])
    reproduction = dict(
        overlap=len(top_deployed & deployed_calls),
        expected=cfg["deployed_calls"],
        exact=top_deployed == deployed_calls,
    )
    print(json.dumps(dict(l4_reproduction=reproduction)), flush=True)

    # The fitted confusion matrix's point MCC is not the published bootstrap mean, so
    # deltas are taken against the point value and the published anchor is reported too.
    point_mcc = live_mcc(state, 0, 0.0)
    anchor = state["mcc"] - point_mcc
    candidates, prediction_rows = [], []
    sources = {
        "l4_rank_mean": ordered["molecule_name"].to_list(),
        "pool_lgbm": combined.sort(["pool_score", "compound_id"], descending=[True, False])[
            "molecule_name"
        ].to_list(),
    }
    for name, counts in cfg["candidates"].items():
        for count in counts:
            calls = set(sources[name][:count])
            added = calls - deployed_calls
            removed = deployed_calls - calls
            scored_share = state["scored_rows"] / 750
            added_scored = len(added) * scored_share
            entry = dict(
                candidate=f"{name}_top{count}",
                score_column=name,
                calls=count,
                added=len(added),
                removed=len(removed),
                changed=len(added) + len(removed),
                monotone_extension=len(removed) == 0,
                expected_added_scored=added_scored,
                readable=len(added) + len(removed) >= cfg["readable_changed_calls"],
            )
            candidates.append(entry)
            if entry["monotone_extension"]:
                for precision in cfg["marginal_precision_scenarios"]:
                    value = live_mcc(state, added_scored, precision)
                    prediction_rows.append(
                        dict(
                            candidate=entry["candidate"],
                            marginal_precision=precision,
                            predicted_point_mcc=value,
                            predicted_anchored_mcc=value + anchor,
                            predicted_delta=value - point_mcc,
                            expected_added_true_positives=added_scored * precision,
                        )
                    )
    # Offline point estimate of the marginal precision of each added rank band, from the
    # frozen held-fold predictions of the same two score columns.
    oof = pl.read_parquet(source / inputs["testlike_oof"])
    band_rows = []
    for split in oof["split_id"].unique().to_list():
        ranked_oof = (
            oof.filter((pl.col("split_id") == split) & pl.col("model").is_in(cfg["l4_components"]))
            .with_columns(
                (
                    (pl.col("score").rank("ordinal").over(["fold", "model"]) - 0.5)
                    / pl.len().over(["fold", "model"])
                ).alias("pct")
            )
            .group_by(["fold", "compound_id"])
            .agg(pl.col("pct").mean().alias("l4_score"))
        )
        roles = oof.filter(
            (pl.col("split_id") == split) & (pl.col("model") == cfg["l4_components"][0])
        ).select("compound_id", "fold", "pool_role", "y_official")
        frame = ranked_oof.join(roles, on=["compound_id", "fold"], validate="1:1")
        frame = frame.with_columns(
            (
                pl.col("l4_score").rank("ordinal", descending=True).over("fold")
                / pl.len().over("fold")
            ).alias("position")
        )
        counts = sorted({0, *cfg["candidates"]["l4_rank_mean"]})
        # Mixture weights are per-row stratum weights; their shared scale cancels in a ratio.
        phi = cfg["identified_phi"]
        weight_official = phi / int((frame["pool_role"] == "official").sum())
        weight_pseudo = (1 - phi) / int((frame["pool_role"] == "pseudo_negative").sum())
        for lower, upper in zip(counts, counts[1:], strict=False):
            band = frame.filter(
                (pl.col("position") > lower / 750) & (pl.col("position") <= upper / 750)
            )
            official = band.filter(pl.col("pool_role") == "official")
            pseudo = band.filter(pl.col("pool_role") == "pseudo_negative")
            positives = float(official["y_official"].sum())
            scored = official.height + pseudo.height
            denominator = weight_official * official.height + weight_pseudo * pseudo.height
            mixture = weight_official * positives / denominator if denominator > 0 else None
            band_rows.append(
                dict(
                    split_id=split,
                    band=f"{lower + 1}-{upper}",
                    live_break_even=cfg["live_break_even_precision"],
                    rows=band.height,
                    official_rows=official.height,
                    official_positives=int(positives),
                    official_precision=positives / official.height if official.height else None,
                    pseudo_rows=pseudo.height,
                    scored_rows=scored,
                    mixture_precision=mixture,
                )
            )
    pl.DataFrame(band_rows).write_csv(out / "marginal_precision.csv")
    pl.DataFrame(candidates).write_csv(out / "candidate_vectors.csv")
    pl.DataFrame(prediction_rows).write_csv(out / "predictions.csv")

    # Resolution: the public metric moves in steps of roughly one scored true positive.
    unit = dict(
        one_true_positive=live_mcc(state, 0, 0.0)
        - mcc_from_counts(
            state["true_positives"] - 1,
            state["scored_calls"] - state["true_positives"] + 1,
            state["scored_positives"] - state["true_positives"] + 1,
            state["scored_rows"]
            - state["scored_calls"]
            - state["scored_positives"]
            + state["true_positives"]
            - 1,
        ),
        one_false_positive=mcc_from_counts(
            state["true_positives"],
            state["scored_calls"] - state["true_positives"] + 1,
            state["scored_positives"] - state["true_positives"],
            state["scored_rows"]
            - state["scored_calls"]
            - 1
            - state["scored_positives"]
            + state["true_positives"],
        )
        - state["mcc"],
    )
    summary = dict(
        config=cfg,
        endpoint=ENDPOINT,
        l4_reproduction=reproduction,
        marginal_precision=band_rows,
        candidates=candidates,
        predictions=prediction_rows,
        metric_units=unit,
        point_mcc=point_mcc,
        published_minus_point_mcc=anchor,
        preregistration=cfg["preregistration"],
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "The Live state is a reconstruction from published aggregates, not disclosed labels",
            "Added ranks are assumed to be scored at the overall scored fraction",
            "No candidate, prepared submission or upload is created by this run",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"], unit=unit)), flush=True)


if __name__ == "__main__":
    main()
