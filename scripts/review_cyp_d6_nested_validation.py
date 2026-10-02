"""Audit completed nested-policy validation, including unknown-prior sensitivity."""

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
    frame = pl.read_parquet(inp / "nested-oof.parquet")
    source = pl.read_parquet(inp / "source-audit.parquet")
    joined = frame.join(
        source.select("compound_id", "regime", "training_weight"), on="compound_id", how="left"
    )
    assert joined.height == frame.height
    summaries, strata = [], []
    for keys, sub in joined.group_by("split_id", "seed", "model", "policy", maintain_order=True):
        key = dict(zip(["split_id", "seed", "model", "policy"], keys, strict=True))
        observed = sub.filter(pl.col("label").is_finite())
        y, call = observed["label"].to_numpy(), observed["call"].to_numpy()
        summaries.append(
            {
                **key,
                "labelled_mcc": float(matthews_corrcoef(y, call)),
                "full_pool_rate": float(sub["call"].mean()),
                "true_positive_rate": float(call[y == 1].mean()),
                "false_positive_rate": float(call[y == 0].mean()),
                **{f"projected_prior_{p:.2f}": projected(y, call, p) for p in cfg["priors"]},
            }
        )
        for column in ["regime", "training_weight"]:
            for value, segment in observed.group_by(column, maintain_order=True):
                truth, calls = segment["label"].to_numpy(), segment["call"].to_numpy()
                strata.append(
                    {
                        **key,
                        "stratum_type": column,
                        "stratum": str(value[0]),
                        "n": len(segment),
                        "positives": int(truth.sum()),
                        "calls": int(calls.sum()),
                        "tp": int(((truth == 1) & calls).sum()),
                        "fp": int(((truth == 0) & calls).sum()),
                        "mcc": float(matthews_corrcoef(truth, calls))
                        if len(np.unique(truth)) > 1
                        else None,
                    }
                )
    pl.DataFrame(summaries).write_csv(out / "actual-policy-prior-sensitivity.csv")
    pl.DataFrame(strata).write_csv(out / "label-regime-decision-audit.csv")
    # Same molecule's call agreement across repeated training seeds.
    agreement = []
    for keys, sub in frame.group_by("split_id", "model", "policy", maintain_order=True):
        a = sub.filter(pl.col("seed") == cfg["seeds"][0]).sort("compound_id")
        b = sub.filter(pl.col("seed") == cfg["seeds"][1]).sort("compound_id")
        assert a["compound_id"].equals(b["compound_id"])
        ca, cb = a["call"].to_numpy(), b["call"].to_numpy()
        union = (ca | cb).sum()
        agreement.append(
            dict(zip(["split_id", "model", "policy"], keys, strict=True))
            | dict(
                changed=int((ca != cb).sum()),
                n=len(a),
                positive_jaccard=float((ca & cb).sum() / union) if union else 1.0,
                score_correlation=float(np.corrcoef(a["score"], b["score"])[0, 1]),
            )
        )
    pl.DataFrame(agreement).write_csv(out / "seed-call-agreement.csv")
    cuts = pl.read_csv(inp / "threshold-scale-audit.csv")
    quantiles = cuts.group_by("split_id", "model").agg(
        pl.col("threshold").min().alias("threshold_min"),
        pl.col("threshold").median().alias("threshold_median"),
        pl.col("threshold").max().alias("threshold_max"),
        (pl.col("refit_valid_q90") - pl.col("inner_score_q90"))
        .median()
        .alias("median_q90_scale_shift"),
        (pl.col("refit_valid_q50") - pl.col("inner_score_q50"))
        .median()
        .alias("median_q50_scale_shift"),
    )
    quantiles.write_csv(out / "threshold-scale-summary.csv")
    audit = pl.read_csv(inp / "leakage-and-operating-points.csv")
    assert not audit["outer_test_used_in_fit_or_cut"].any()
    if (inp / "comparison-oof.parquet").exists():
        from validate_cyp_d6_bag_policy import paired_intervals

        reference = pl.read_parquet(inp / "comparison-oof.parquet").filter(
            pl.col("model") == "fresh_d6_network"
        )
        current = frame.filter(pl.col("model") == "fresh_d6_network").with_columns(
            pl.lit("fresh_d6_network_20epochs").alias("model")
        )
        paired_intervals(
            pl.concat([reference, current]),
            out,
            {
                "seeds": cfg["seeds"],
                "bootstrap_samples": 1000,
                "split_strategies": ["umap", "scaffold"],
                "models": ["fresh_d6_network", "fresh_d6_network_20epochs"],
            },
        )
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "fit_boundary_checks_pass": True,
                "interpretation": "Prior scenarios are hypothetical; supplied-label regime "
                "slices do not reveal the hidden scored mask. Same outer folds have been "
                "inspected previously, so this is conditional protocol validation, "
                "not new independent model-family selection evidence.",
                "no_upload": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
