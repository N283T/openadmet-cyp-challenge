"""Audit recommendation evidence and decision-policy confounding using saved predictions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import calls_at_rate
from sklearn.metrics import matthews_corrcoef


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    data = pl.read_parquet(inp / "combined-oof.parquet")
    previous = pl.read_csv(inp / "matched-rank-fractions.csv")
    history = pl.read_csv(inp / "historical-concordance.csv")
    history.write_csv(out / "prior-historical-concordance.csv")
    rows = []
    for split in ["umap", "scaffold"]:
        parts = [
            data.filter((pl.col("split_id") == split) & (pl.col("model") == model)).sort(
                "compound_id"
            )
            for model in [cfg["legacy_model"], cfg["selected_model"]]
        ]
        old, new = parts
        assert old.height == new.height == 6146
        assert old.select("compound_id", "label", "fold", "group").equals(
            new.select("compound_id", "label", "fold", "group")
        )
        known = np.isfinite(old["label"].to_numpy())
        truth = old["label"].to_numpy()[known]
        old_trained = old["call"].to_numpy()[known]
        old_rank = calls_at_rate(
            old["score"].to_numpy(), old["fold"].to_numpy(), cfg["historical_fraction"]
        )[known]
        new_trained = new["call"].to_numpy()[known]
        rank_reference = previous.filter(
            (pl.col("split_id") == split)
            & (pl.col("model") == cfg["legacy_model"])
            & ((pl.col("requested_fraction") - cfg["historical_fraction"]).abs() < 1e-12)
        )
        assert rank_reference.height == 1
        np.testing.assert_allclose(
            matthews_corrcoef(truth, old_rank), rank_reference["mcc"][0], atol=1e-12
        )
        groups = old["group"].to_numpy()[known]
        members = [np.flatnonzero(groups == g) for g in sorted(set(groups))]
        for reference_name, reference in [
            ("legacy_new_trained_threshold", old_trained),
            ("legacy_historical_rank_fraction", old_rank),
        ]:
            rng = np.random.default_rng(20261002)
            deltas = []
            for _ in range(cfg["bootstrap_samples"]):
                indices = np.concatenate(
                    [members[i] for i in rng.integers(0, len(members), len(members))]
                )
                deltas.append(
                    matthews_corrcoef(truth[indices], new_trained[indices])
                    - matthews_corrcoef(truth[indices], reference[indices])
                )
            reference_mcc = float(matthews_corrcoef(truth, reference))
            selected_mcc = float(matthews_corrcoef(truth, new_trained))
            rows.append(
                dict(
                    split_id=split,
                    comparator=reference_name,
                    known_n=int(known.sum()),
                    reference_mcc=reference_mcc,
                    selected_mcc=selected_mcc,
                    delta_mcc=selected_mcc - reference_mcc,
                    ci_low=float(np.quantile(deltas, 0.025)),
                    ci_high=float(np.quantile(deltas, 0.975)),
                    changed_known_calls=int(np.sum(reference != new_trained)),
                )
            )
    pl.DataFrame(rows).write_csv(out / "decision-reference-comparisons.csv")
    full_pool = history.filter(pl.col("scope").str.starts_with("complete-pool"))
    assert full_pool.height == 4 and not full_pool["concordant"].any()
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                comparisons=rows,
                prior_full_pool_concordance=full_pool.to_dicts(),
                conclusion="Whole-recipe same-threshold comparison was not validation "
                "of prospective Live model selection. Previously measured historical "
                "discordance remained unresolved. UMAP advantage disappears against "
                "the legacy historical-fraction diagnostic; scaffold advantage remains. "
                "This isolates a comparator-policy weakness, not the exact molecular "
                "cause of the Live failure or a validated replacement selector.",
                no_training=True,
                no_new_candidate=True,
                no_upload=True,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
