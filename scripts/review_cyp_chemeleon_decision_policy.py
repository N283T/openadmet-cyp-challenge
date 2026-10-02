"""Apply the predeclared decision-policy choice to exact two-seed averaged outer predictions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import projected
from sklearn.metrics import matthews_corrcoef


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    pred = pl.read_parquet(inp / "predictions.parquet")
    original = pl.read_parquet(inp / "original-predictions.parquet").filter(
        (pl.col("model") == "chemeleon_primary_network")
        & (pl.col("policy") == "centered_fold_mean")
    )
    current = pred.filter(pl.col("policy") == "centered_fold_mean")
    keys = ["split_id", "seed", "fold", "compound_id"]
    original, current = original.sort(keys), current.sort(keys)
    assert original.select(keys).equals(current.select(keys))
    np.testing.assert_array_equal(original["score"].to_numpy(), current["score"].to_numpy())
    cuts = pl.read_csv(inp / "cuts.csv")
    frames = []
    for policy, column in [
        ("centered_fold_mean", "threshold"),
        ("prior_mean_centered", "prior_mean_cut"),
    ]:
        threshold = (
            cuts.group_by("split_id", "seed", "outer_fold", "model")
            .agg(pl.col(column).mean().alias("cut"))
            .rename({"outer_fold": "fold"})
        )
        sub = (
            pred.filter(pl.col("policy") == policy)
            .join(threshold, on=["split_id", "seed", "fold", "model"], validate="m:1")
            .with_columns((pl.col("score") - pl.col("cut")).alias("margin"))
        )
        assert sub.filter((pl.col("margin") > 0) != pl.col("call")).height == 0
        averaged = (
            sub.group_by("split_id", "compound_id", "fold", "group", "endpoint")
            .agg(
                pl.col("label").first(),
                pl.col("score").mean(),
                pl.col("margin").mean(),
                pl.col("seed").n_unique().alias("seeds"),
            )
            .with_columns(pl.lit(policy).alias("policy"), (pl.col("margin") > 0).alias("call"))
        )
        assert (averaged["seeds"] == 2).all()
        frames.append(averaged)
    combined = pl.concat(frames)
    combined.write_parquet(out / "combined-policy-predictions.parquet")
    rows, decision = [], {}
    for split in ["umap", "scaffold"]:
        observed = combined.filter((pl.col("split_id") == split) & pl.col("label").is_finite())
        base = observed.filter(pl.col("policy") == "centered_fold_mean").sort("compound_id")
        other = observed.filter(pl.col("policy") == "prior_mean_centered").sort("compound_id")
        assert base["compound_id"].equals(other["compound_id"])
        y, a, b = (base["label"].to_numpy(), base["call"].to_numpy(), other["call"].to_numpy())
        deltas = []
        for prior in [None, 0.03, 0.07, 0.12]:
            before = matthews_corrcoef(y, a) if prior is None else projected(y, a, prior)
            after = matthews_corrcoef(y, b) if prior is None else projected(y, b, prior)
            rows.append(
                dict(
                    split_id=split,
                    assumed_prior=prior,
                    observed_cut_mcc=float(before),
                    mean_prior_cut_mcc=float(after),
                    delta=float(after - before),
                )
            )
            if prior is not None:
                deltas.append(after - before)
        decision[split] = bool(deltas[1] > 0 and np.mean(deltas) > 0)
    selected = "prior_mean_centered" if all(decision.values()) else "centered_fold_mean"
    pl.DataFrame(rows).write_csv(out / "combined-policy-metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                selected_policy=selected,
                selection_rule_passed=decision,
                metrics=rows,
                original_scores_exactly_reproduced=True,
                individual_calls_exactly_reconstructed=True,
                no_upload=True,
                limits="Predeclared decision on conditional class-label-shift scenarios, "
                "not known blind prevalence or identified weak-primary outcomes. "
                "No candidate count was used in the choice.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
