"""Compare exact averaged fresh policies and frozen historical full-pool controls."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from report_cyp_primary_selection_controls import paired_delta
from review_cyp_d6_causal_probe import calls_at_rate, projected
from sklearn.metrics import average_precision_score, matthews_corrcoef


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    modern = pl.read_parquet(inp / "combined-oof.parquet")
    curves, contrasts, stability = [], [], []
    for split in ["umap", "scaffold"]:
        data = modern.filter(pl.col("split_id") == split)
        for model in ["fresh_d6_network", "chemeleon_class_network"]:
            a = (
                data.filter(pl.col("model") == "chemeleon_primary_network")
                .filter(pl.col("label").is_finite())
                .with_columns(pl.lit(0).alias("seed"))
            )
            b = (
                data.filter(pl.col("model") == model)
                .filter(pl.col("label").is_finite())
                .with_columns(pl.lit(0).alias("seed"))
            )
            contrasts.append(
                dict(
                    split_id=split,
                    model="chemeleon_primary_network",
                    comparator=model,
                    **paired_delta(a, b, cfg["bootstrap_samples"]),
                )
            )
        historical = (
            pl.read_parquet(inp / f"control-{split}.parquet")
            .filter(pl.col("model") == "sep15_control_raw_oof_analog")
            .rename({"y_true": "label"})
        )
        baseline = historical.sort("compound_id")
        frames = [(m, data.filter(pl.col("model") == m)) for m in data["model"].unique().sort()]
        frames.append(("historical_control_analog", historical))
        for name, frame in frames:
            frame = frame.sort("compound_id")
            assert frame["compound_id"].equals(baseline["compound_id"])
            assert frame["fold"].equals(baseline["fold"])
            np.testing.assert_allclose(frame["label"], baseline["label"], equal_nan=True)
            observed = np.isfinite(frame["label"].to_numpy())
            y = frame["label"].to_numpy()[observed]
            for rate in cfg["rates"]:
                calls = calls_at_rate(frame["score"].to_numpy(), frame["fold"].to_numpy(), rate)
                for prior in cfg["priors"]:
                    curves.append(
                        dict(
                            split_id=split,
                            model=name,
                            full_pool_rank_rate=rate,
                            assumed_prior=prior,
                            projected_mcc=projected(y, calls[observed], prior),
                            mcc=float(matthews_corrcoef(y, calls[observed])),
                            ap=float(
                                average_precision_score(y, frame["score"].to_numpy()[observed])
                            ),
                            labelled_calls=int(calls[observed].sum()),
                            unlabelled_calls=int(calls[~observed].sum()),
                        )
                    )
    blind = pl.read_parquet(inp / "blind-folds.parquet")
    decisions = []
    for (seed,), sub in blind.group_by("seed"):
        result = sub.group_by("compound_id").agg(pl.col("margin").mean()).sort("compound_id")
        decisions.append((seed, result))
    a, b = sorted(decisions, key=lambda x: x[0])
    assert a[1]["compound_id"].equals(b[1]["compound_id"])
    ca, cb = (r["margin"].to_numpy() > 0 for _, r in [a, b])
    stability.append(
        dict(
            seed_a=a[0],
            seed_b=b[0],
            a_positives=int(ca.sum()),
            b_positives=int(cb.sum()),
            changed=int((ca != cb).sum()),
            positive_jaccard=float((ca & cb).sum() / max((ca | cb).sum(), 1)),
        )
    )
    pl.DataFrame(curves).write_csv(out / "historical-full-pool-curves.csv")
    pl.DataFrame(contrasts).write_csv(out / "paired-fresh-policy-contrasts.csv")
    pl.DataFrame(stability).write_csv(out / "blind-seed-agreement.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                contrasts=contrasts,
                blind_seed_agreement=stability,
                no_upload=True,
                limits="Historical recipe is a frozen analog, "
                "not exact nested historical model selection. "
                "Rate curves are sensitivity diagnostics, not a fitted call budget. "
                "Conditional group intervals omit protocol and training uncertainty. "
                "Missing weak-primary outcomes and unknown blind scored mask remain unidentified.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
