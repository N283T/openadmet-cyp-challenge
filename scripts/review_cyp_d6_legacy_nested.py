"""Compare complete legacy and fresh candidates on identical untouched outer records."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import projected
from sklearn.metrics import average_precision_score, matthews_corrcoef


def combine(frame, expected_seeds):
    frame = frame.filter(pl.col("policy") == "centered_fold_mean")
    if "margin" not in frame.columns:
        raise ValueError("A trained-cut margin is required")
    assert frame.filter((pl.col("margin") > 0) != pl.col("call")).height == 0
    result = frame.group_by("split_id", "model", "compound_id", "fold", "group", "endpoint").agg(
        pl.col("label").first(),
        pl.col("label").n_unique().alias("label_versions"),
        pl.col("score").mean(),
        pl.col("margin").mean(),
        pl.col("seed").n_unique().alias("seeds"),
    )
    assert (result["seeds"] == expected_seeds).all()
    assert (result["label_versions"] == 1).all()
    return result.with_columns((pl.col("margin") > 0).alias("call"))


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    assert json.loads((inp / "legacy-result.json").read_text())["exit_code"] == 0
    legacy = pl.read_parquet(inp / "legacy.parquet")
    fresh = pl.read_parquet(inp / "fresh.parquet").filter(pl.col("policy") == "centered_fold_mean")
    cuts = (
        pl.read_csv(inp / "fresh-cuts.csv")
        .group_by("split_id", "seed", "outer_fold", "model")
        .agg(pl.col("threshold").mean().alias("cut"))
        .rename({"outer_fold": "fold"})
    )
    fresh = fresh.join(cuts, on=["split_id", "seed", "fold", "model"], validate="m:1")
    fresh = fresh.with_columns((pl.col("score") - pl.col("cut")).alias("margin"))
    combined = pl.concat([combine(legacy, 2), combine(fresh, 2)], how="vertical_relaxed")
    combined.write_parquet(out / "combined-oof.parquet")
    rows, contrasts, rate_rows = [], [], []
    for (split,), block in combined.group_by("split_id"):
        frames = {model: frame.sort("compound_id") for (model,), frame in block.group_by("model")}
        old = frames[cfg["legacy_model"]]
        assert old.height == 6146 and old["compound_id"].n_unique() == 6146
        keys = ["compound_id", "endpoint", "split_id", "fold", "group", "label"]
        for model, frame in frames.items():
            assert frame.select(keys).equals(old.select(keys))
            observed = frame.filter(pl.col("label").is_finite())
            truth, calls = observed["label"].to_numpy(), observed["call"].to_numpy()
            rows.append(
                dict(
                    split_id=split,
                    model=model,
                    ap=float(average_precision_score(truth, observed["score"])),
                    mcc=float(matthews_corrcoef(truth, calls)),
                    full_pool_rate=float(frame["call"].mean()),
                    labelled_rate=float(observed["call"].mean()),
                    **{
                        f"projected_mcc_prior_{p}": projected(truth, calls, p)
                        for p in [0.03, 0.07, 0.12]
                    },
                )
            )
            # Diagnostic matched rank fractions are not trained decisions or budgets.
            for rate in cfg["diagnostic_rank_fractions"]:
                ranked_calls = np.zeros(frame.height, dtype=bool)
                membership = frame["fold"].to_numpy()
                values = frame["score"].to_numpy()
                for fold in range(5):
                    indices = np.flatnonzero(membership == fold)
                    count = max(1, int(round(rate * len(indices))))
                    order = np.argsort(-values[indices], kind="stable")
                    ranked_calls[indices[order[:count]]] = True
                labelled = np.isfinite(frame["label"].to_numpy())
                rate_rows.append(
                    dict(
                        split_id=split,
                        model=model,
                        requested_fraction=rate,
                        actual_fraction=float(ranked_calls.mean()),
                        mcc=float(matthews_corrcoef(truth, ranked_calls[labelled])),
                        projected_mcc_prior_007=projected(truth, ranked_calls[labelled], 0.07),
                    )
                )
            if model == cfg["legacy_model"]:
                continue
            prior_old = old.filter(pl.col("label").is_finite())
            a, b = calls, prior_old["call"].to_numpy()
            groups = observed["group"].to_numpy()
            members = [np.flatnonzero(groups == g) for g in sorted(set(groups))]
            rng = np.random.default_rng(20261002)
            deltas = []
            for _ in range(cfg["bootstrap_samples"]):
                indices = np.concatenate(
                    [members[i] for i in rng.integers(0, len(members), len(members))]
                )
                deltas.append(
                    matthews_corrcoef(truth[indices], a[indices])
                    - matthews_corrcoef(truth[indices], b[indices])
                )
            contrasts.append(
                dict(
                    split_id=split,
                    model=model,
                    comparator=cfg["legacy_model"],
                    delta_mcc=float(matthews_corrcoef(truth, a) - matthews_corrcoef(truth, b)),
                    ci_low=float(np.quantile(deltas, 0.025)),
                    ci_high=float(np.quantile(deltas, 0.975)),
                    changed_known_calls=int(np.sum(a != b)),
                )
            )
    pl.DataFrame(rows).write_csv(out / "policy-comparison.csv")
    pl.DataFrame(contrasts).write_csv(out / "paired-contrasts.csv")
    pl.DataFrame(rate_rows).write_csv(out / "matched-rank-fractions.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                policies=rows,
                paired_contrasts=contrasts,
                exact_identity_fold_label_match=True,
                limits="Conditional intervals omit fitting and recipe selection uncertainty. "
                "This compares held-out decisions on observed D6 labels, not missing "
                "weak-primary outcomes or prospective leaderboard ordering.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
