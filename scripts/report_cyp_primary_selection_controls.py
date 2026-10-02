"""Measure selection reversal on known held-group outcomes and compare gate mechanisms."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import calls_at_rate
from sklearn.metrics import average_precision_score, matthews_corrcoef


def paired_delta(a, b, samples):
    joined = a.select("compound_id", "seed", "group", "label", "call").join(
        b.select("compound_id", "seed", pl.col("call").alias("other")),
        on=["compound_id", "seed"],
        validate="1:1",
    )
    assert joined.height == a.height == b.height
    y, calls, other, seeds = (joined[c].to_numpy() for c in ["label", "call", "other", "seed"])
    groups = joined["group"].to_numpy()
    group_ids = sorted(set(groups))
    members = [np.flatnonzero(groups == group) for group in group_ids]
    seed_ids = sorted(set(seeds))

    def difference(idx):
        return np.mean(
            [
                matthews_corrcoef(y[ix], calls[ix]) - matthews_corrcoef(y[ix], other[ix])
                for seed in seed_ids
                if len(ix := idx[seeds[idx] == seed])
            ]
        )

    rng = np.random.default_rng(20261002)
    draws = [
        difference(
            np.concatenate([members[i] for i in rng.integers(0, len(members), len(members))])
        )
        for _ in range(samples)
    ]
    return dict(
        delta=float(difference(np.arange(len(y)))),
        ci_low=float(np.quantile(draws, 0.025)),
        ci_high=float(np.quantile(draws, 0.975)),
    )


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    summaries, contrasts, stability, coverage, rate_metrics = [], [], [], [], []
    for case in cfg["cases"]:
        frame = pl.read_parquet(inp / case["predictions"]).filter(
            pl.col("policy") == "centered_fold_mean"
        )
        source = pl.read_parquet(inp / case["source"])
        labelled = source.filter(pl.col("label").is_finite())
        for (split, seed, model), sub in frame.group_by("split_id", "seed", "model"):
            sub = sub.sort("compound_id")
            observed = sub["label"].is_finite().to_numpy()
            for rate in [0.03, 0.07, 0.12, 0.33]:
                calls = calls_at_rate(sub["score"].to_numpy(), sub["fold"].to_numpy(), rate)
                truth = sub["label"].to_numpy()[observed]
                rate_metrics.append(
                    dict(
                        case=case["name"],
                        split_id=split,
                        seed=seed,
                        model=model,
                        full_pool_rank_rate=rate,
                        mcc=float(matthews_corrcoef(truth, calls[observed])),
                        labelled_calls=int(calls[observed].sum()),
                        labelled_true_positive=int((calls[observed] & (truth == 1)).sum()),
                    )
                )
        if "selected" in source.columns:
            for selected, sub in labelled.group_by("selected"):
                coverage.append(
                    dict(
                        case=case["name"],
                        selected=selected[0],
                        n=sub.height,
                        positives=int(sub["label"].sum()),
                    )
                )
        for scope in ["all_known", "selected", "weak_or_missing"]:
            if scope != "all_known" and not case["selection_control"]:
                continue
            data = frame.filter(pl.col("label").is_finite())
            if scope != "all_known":
                data = data.filter(pl.col("selected") == (scope == "selected"))
            if not data.height:
                continue
            for (split, seed, model), sub in data.group_by("split_id", "seed", "model"):
                y, score, calls = (sub[c].to_numpy() for c in ["label", "score", "call"])
                summaries.append(
                    dict(
                        case=case["name"],
                        scope=scope,
                        split_id=split,
                        seed=seed,
                        model=model,
                        n=len(y),
                        positives=int(y.sum()),
                        ap=float(average_precision_score(y, score)),
                        mcc=float(matthews_corrcoef(y, calls)),
                        call_rate=float(calls.mean()),
                        false_positive=int(((y == 0) & calls).sum()),
                        false_negative=int(((y == 1) & ~calls).sum()),
                    )
                )
            for split in data["split_id"].unique().sort():
                for a, b in case["contrasts"]:
                    x = data.filter((pl.col("split_id") == split) & (pl.col("model") == a))
                    z = data.filter((pl.col("split_id") == split) & (pl.col("model") == b))
                    contrasts.append(
                        dict(
                            case=case["name"],
                            scope=scope,
                            split_id=split,
                            model=a,
                            comparator=b,
                            **paired_delta(x, z, cfg["bootstrap_samples"]),
                        )
                    )
        for (split, model), sub in frame.group_by("split_id", "model"):
            seeds = sub["seed"].unique().sort().to_list()
            a = sub.filter(pl.col("seed") == seeds[0]).sort("compound_id")
            b = sub.filter(pl.col("seed") == seeds[1]).sort("compound_id")
            assert a["compound_id"].equals(b["compound_id"])
            ac, bc = a["call"].to_numpy(), b["call"].to_numpy()
            stability.append(
                dict(
                    case=case["name"],
                    split_id=split,
                    model=model,
                    changed=int((ac != bc).sum()),
                    positive_jaccard=float((ac & bc).sum() / max((ac | bc).sum(), 1)),
                    score_correlation=float(np.corrcoef(a["score"], b["score"])[0, 1]),
                )
            )
    pl.DataFrame(summaries).write_csv(out / "population-metrics.csv")
    pl.DataFrame(contrasts).write_csv(out / "paired-population-contrasts.csv")
    pl.DataFrame(stability).write_csv(out / "seed-stability.csv")
    pl.DataFrame(coverage).write_csv(out / "known-label-coverage.csv")
    pl.DataFrame(rate_metrics).write_csv(out / "matched-full-pool-rate-metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                contrasts=contrasts,
                no_upload=True,
                unknown_classes_never_imputed=True,
                limits="Conditional group resampling of saved fits; "
                "protocol and training uncertainty omitted. "
                "A4 control establishes a possible selection mechanism, "
                "not the realized cause of D6 Live failures.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
