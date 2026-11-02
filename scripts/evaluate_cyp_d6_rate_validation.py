"""Paired validation of the D6 call count for the deployed L4 rank-mean score.

The blind D6 set is an analog expansion around 1A2/2C9/3A4 hits (D6 was not an anchor), is dense
across enzymes while the training D6 labels are sparse, and is scored only where a label can be
assigned. Three populations are therefore used, each as a paired grouped bootstrap of
MCC(k calls) - MCC(49 calls) on the frozen held-fold L4 scores:

* Metric O: officially labelled rows, classes reweighted to a target scored prevalence;
* Metric M: mixture of the labelled stratum (mass phi) and the pseudo-negative stratum;
* analog pool: the 608 pseudo analog-expansion molecules ranked together, labelled rows only.

Nothing is fitted, no candidate is created and no API is contacted.
"""

import json
import os
from pathlib import Path

import evaluate_cyp_d6_official_metric as om
import numpy as np
import polars as pl


def l4_frame(oof, components):
    ranked = (
        oof.filter(pl.col("model").is_in(components))
        .with_columns(
            (
                (pl.col("score").rank("ordinal").over(["split_id", "fold", "model"]) - 0.5)
                / pl.len().over(["split_id", "fold", "model"])
            ).alias("pct")
        )
        .group_by(["split_id", "fold", "compound_id"])
        .agg(pl.col("pct").mean().alias("score"), pl.len().alias("members"))
    )
    assert ranked["members"].min() == len(components)
    template = oof.filter(pl.col("model") == om.CONTROL).drop("score")
    frame = template.join(
        ranked.drop("members"), on=["split_id", "fold", "compound_id"], how="left"
    )
    assert frame["score"].null_count() == 0
    return frame


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    oof = pl.read_parquet(inp / "testlike-pool-oof.parquet")
    assignments = pl.read_parquet(inp / "assignments.parquet")
    analog = set(pl.read_parquet(inp / "analog-pairs.parquet")["analog_id"].to_list())
    base_rate = cfg["reference_calls"] / 750
    rows = []
    for split in sorted(oof["split_id"].unique().to_list()):
        frame = (
            l4_frame(oof.filter(pl.col("split_id") == split), cfg["l4_components"])
            .sort("compound_id")
            .join(
                assignments.filter(pl.col("strategy") == split).select("compound_id", "group"),
                on="compound_id",
                how="left",
            )
            .with_columns(pl.col(["y_official", "y_testlike"]).fill_null(np.nan))
        )
        assert frame["group"].null_count() == 0
        score, fold = frame["score"].to_numpy(), frame["fold"].to_numpy()
        groups = frame["group"].to_numpy().astype(str)
        role = frame["pool_role"].to_numpy()
        official = role == "official"
        pseudo = role == "pseudo_negative"
        y_official, y_testlike = frame["y_official"].to_numpy(), frame["y_testlike"].to_numpy()
        in_pool = np.array([i in analog for i in frame["compound_id"].to_list()])

        populations = []
        for prior in cfg["priors"]:
            populations.append(
                (
                    "O",
                    prior,
                    "dev",
                    om.official_weights(y_official, official, prior),
                    y_official,
                    None,
                )
            )
        for phi in cfg["phis"]:
            populations.append(
                ("M", phi, "dev", om.mixture_weights(official, pseudo, phi), y_testlike, None)
            )
        for prior in cfg["priors"]:
            populations.append(
                (
                    "O",
                    prior,
                    "analog_pool",
                    om.official_weights(y_official, official & in_pool, prior),
                    y_official,
                    in_pool,
                )
            )
        for metric, parameter, population, weights, labels, pool in populations:

            def calls_for(count, pool=pool, score=score, fold=fold):
                rate = count / 750
                if pool is None:
                    return om.top_calls(score, fold, rate)
                calls = np.zeros(len(score), dtype=bool)
                members = np.flatnonzero(pool)
                top = members[
                    np.argsort(-score[members], kind="stable")[: round(rate * len(members))]
                ]
                calls[top] = True
                return calls

            reference = calls_for(cfg["reference_calls"])
            point_reference = om.weighted_mcc(labels, reference, weights)
            for count in cfg["counts"]:
                calls = calls_for(count)
                result = om.group_bootstrap(
                    labels, groups, weights, calls, reference, cfg["bootstrap_samples"], cfg["seed"]
                )
                rows.append(
                    dict(
                        split=split,
                        population=population,
                        metric=metric,
                        parameter=parameter,
                        calls=count,
                        mcc_reference=point_reference,
                        mcc_calls=om.weighted_mcc(labels, calls, weights),
                        **result,
                    )
                )
    table = pl.DataFrame(rows)
    table.write_csv(out / "rate-validation.csv")
    summary = table.filter(pl.col("calls") == cfg["decision_calls"]).select(
        "split",
        "population",
        "metric",
        "parameter",
        "delta",
        "ci_low",
        "ci_high",
        "p_delta_gt_0",
        "molecules",
    )
    print(summary)
    (out / "summary.json").write_text(json.dumps(summary.to_dicts(), indent=2) + "\n")
    del base_rate


if __name__ == "__main__":
    main()
