"""Validate equal-seed margin averaging using complete untouched-outer fit records."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import projected
from sklearn.metrics import average_precision_score, matthews_corrcoef


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    pred = pl.read_parquet(inp / "nested-oof.parquet").filter(
        pl.col("policy") == "centered_fold_mean"
    )
    cuts = (
        pl.read_csv(inp / "cuts.csv")
        .group_by("split_id", "seed", "outer_fold", "model")
        .agg(pl.col("threshold").mean().alias("cut"))
        .rename({"outer_fold": "fold"})
    )
    pred = pred.join(cuts, on=["split_id", "seed", "fold", "model"], validate="m:1").with_columns(
        (pl.col("score") - pl.col("cut")).alias("margin")
    )
    assert pred.filter((pl.col("margin") > 0) != pl.col("call")).height == 0
    assert pred["cut"].is_finite().all()
    combined = (
        pred.group_by("split_id", "model", "compound_id", "fold", "group", "endpoint")
        .agg(
            pl.col("label").first(),
            pl.col("label").n_unique().alias("label_versions"),
            pl.col("score").mean(),
            pl.col("margin").mean(),
            pl.col("seed").n_unique().alias("seeds"),
            pl.col("call").first().alias("seed_a"),
            pl.col("call").last().alias("seed_b"),
        )
        .with_columns((pl.col("margin") > 0).alias("call"))
    )
    assert (combined["seeds"] == 2).all() and (combined["label_versions"] == 1).all()
    combined.write_parquet(out / "combined-oof.parquet")
    rows = []
    for (split, model), frame in combined.group_by("split_id", "model"):
        observed = frame.filter(pl.col("label").is_finite())
        y, calls, a, b, groups = (
            observed[c].to_numpy() for c in ["label", "call", "seed_a", "seed_b", "group"]
        )
        group_ids = sorted(set(groups))
        members = [np.flatnonzero(groups == g) for g in group_ids]
        for prior in [None, 0.03, 0.07, 0.12]:

            def value(y, call, prior=prior):
                return matthews_corrcoef(y, call) if prior is None else projected(y, call, prior)

            rng = np.random.default_rng(20261002)
            deltas = []
            for _ in range(cfg["bootstrap_samples"]):
                idx = np.concatenate(
                    [members[i] for i in rng.integers(0, len(members), len(members))]
                )
                deltas.append(
                    value(y[idx], calls[idx]) - (value(y[idx], a[idx]) + value(y[idx], b[idx])) / 2
                )
            individual_mean = (value(y, a) + value(y, b)) / 2
            rows.append(
                dict(
                    split_id=split,
                    model=model,
                    assumed_prior=prior,
                    combined_mcc=float(value(y, calls)),
                    individual_mcc=float(individual_mean),
                    delta=float(value(y, calls) - individual_mean),
                    ci_low=float(np.quantile(deltas, 0.025)),
                    ci_high=float(np.quantile(deltas, 0.975)),
                    ap=float(average_precision_score(y, observed["score"])),
                    full_pool_rate=float(frame["call"].mean()),
                )
            )
    pl.DataFrame(rows).write_csv(out / "seed-averaging-metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                rows=rows,
                exact_individual_call_reconstruction=True,
                seeds=2,
                no_upload=True,
                limits="Same outer groups, all fits/cuts already exclude outer labels. "
                "Averaging was proposed after observing seed instability. "
                "Conditional intervals omit protocol selection and training uncertainty; "
                "hypothetical priors do not identify absent weak-region outcomes.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
