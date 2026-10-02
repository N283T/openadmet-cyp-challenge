"""Compare genuinely inner-selected prior cuts at matching held-out label-shift scenarios."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import projected


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    data = pl.read_parquet(inp / "nested-oof.parquet")
    rows, intervals = [], []
    for (split, model), frame in data.group_by("split_id", "model"):
        for prior in cfg["priors"]:
            policies = ["centered_fold_mean", f"prior_{prior}_centered"]
            labelled = frame.filter(pl.col("label").is_finite())
            for seed in frame["seed"].unique().sort():
                for policy in policies:
                    all_pool = frame.filter((pl.col("seed") == seed) & (pl.col("policy") == policy))
                    sub = all_pool.filter(pl.col("label").is_finite())
                    rows.append(
                        dict(
                            split_id=split,
                            model=model,
                            assumed_prior=prior,
                            seed=seed,
                            policy=policy,
                            projected_mcc=projected(
                                sub["label"].to_numpy(), sub["call"].to_numpy(), prior
                            ),
                            full_pool_call_rate=float(all_pool["call"].mean()),
                        )
                    )
            joined = (
                labelled.filter(pl.col("policy") == policies[0])
                .select("compound_id", "seed", "group", "label", pl.col("call").alias("baseline"))
                .join(
                    labelled.filter(pl.col("policy") == policies[1]).select(
                        "compound_id", "seed", pl.col("call").alias("prior_call")
                    ),
                    on=["compound_id", "seed"],
                    validate="1:1",
                )
            )
            y, a, b, groups, seeds = (
                joined[c].to_numpy() for c in ["label", "baseline", "prior_call", "group", "seed"]
            )
            group_ids = sorted(set(groups))
            members = [np.flatnonzero(groups == g) for g in group_ids]
            seed_ids = sorted(set(seeds))

            def delta(idx, y=y, a=a, b=b, prior=prior, seed_ids=seed_ids, seeds=seeds):
                return np.mean(
                    [
                        projected(y[ix], b[ix], prior) - projected(y[ix], a[ix], prior)
                        for seed in seed_ids
                        if len(ix := idx[seeds[idx] == seed])
                    ]
                )

            rng = np.random.default_rng(20261002)
            samples = [
                delta(
                    np.concatenate(
                        [members[i] for i in rng.integers(0, len(members), len(members))]
                    )
                )
                for _ in range(cfg["bootstrap_samples"])
            ]
            intervals.append(
                dict(
                    split_id=split,
                    model=model,
                    assumed_prior=prior,
                    delta=float(delta(np.arange(len(y)))),
                    ci_low=float(np.quantile(samples, 0.025)),
                    ci_high=float(np.quantile(samples, 0.975)),
                )
            )
    pl.DataFrame(rows).write_csv(out / "matching-prior-metrics.csv")
    pl.DataFrame(intervals).write_csv(out / "paired-prior-policy-contrasts.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                contrasts=intervals,
                no_upload=True,
                limits="Hypothetical label shift, not known blind prevalence or blind accuracy. "
                "Conditional group intervals exclude fit and protocol uncertainty; "
                "unseen weak-primary "
                "D6 class outcomes remain unidentified.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
