"""Audit positive-call counts for a complete-pool classification OOF model."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.classification_evaluation import fixed_rate_calls, projected_mcc


def load_oof(path: Path, model: str, endpoint: str, split: str, expected_rows: int) -> pl.DataFrame:
    frame = (
        pl.read_parquet(path)
        .filter((pl.col("model") == model) & (pl.col("split_id") == split))
        .sort("compound_id")
    )
    if frame.height != expected_rows or frame["compound_id"].n_unique() != expected_rows:
        raise ValueError(f"Incomplete or duplicated molecular pool: {path}")
    if frame["fold"].n_unique() != 5:
        raise ValueError(f"Expected five outer folds: {path}")
    if not np.isfinite(frame["score"].to_numpy()).all():
        raise ValueError(f"Non-finite scores: {path}")
    if frame["endpoint"].unique().to_list() != [endpoint]:
        raise ValueError(f"Wrong endpoint: {path}")
    return frame


def main() -> None:
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    counts = cfg["counts_per_750"]
    reference = cfg["reference_count"]
    if reference not in counts or len(counts) != len(set(counts)):
        raise ValueError("Reference must occur exactly once in the count grid")
    rates = {count: count / 750 for count in counts}
    rng = np.random.default_rng(cfg["bootstrap_seed"])
    curve_rows = []
    choice_rows = []
    bootstrap_rows = []
    policy_rows = []
    for split, alias in cfg["oof"].items():
        frame = load_oof(
            inputs / alias, cfg["model"], cfg["endpoint"], split, cfg["expected_pool_rows"]
        )
        folds = frame["fold"].to_numpy()
        scores = frame["score"].to_numpy()
        raw_y = frame["y_true"].to_numpy()
        labelled = ~np.isnan(raw_y.astype(float))
        if not labelled.any():
            raise ValueError(f"No scored labels: {split}")
        groups = frame["group"].to_numpy()
        calls = {count: fixed_rate_calls(scores, folds, rate) for count, rate in rates.items()}
        reference_calls = calls[reference]
        unique_groups = np.unique(groups[labelled])
        group_rows = [np.flatnonzero(labelled & (groups == group)) for group in unique_groups]
        samples = [
            np.concatenate(
                [group_rows[index] for index in rng.integers(len(group_rows), size=len(group_rows))]
            )
            for _ in range(cfg["bootstrap_replicates"])
        ]
        for prior in cfg["scored_population_priors"]:
            nested_call = np.zeros(frame.height, dtype=bool)
            for count in counts:
                call = calls[count]
                observed = projected_mcc(raw_y[labelled].astype(int), call[labelled], prior)
                baseline = projected_mcc(
                    raw_y[labelled].astype(int), reference_calls[labelled], prior
                )
                deltas = []
                for rows in samples:
                    y = raw_y[rows].astype(int)
                    delta = projected_mcc(y, call[rows], prior) - projected_mcc(
                        y, reference_calls[rows], prior
                    )
                    if np.isfinite(delta):
                        deltas.append(delta)
                if not deltas:
                    raise ValueError("No valid bootstrap samples")
                curve_rows.append(
                    {
                        "split_id": split,
                        "prior": prior,
                        "count_per_750": count,
                        "rate": rates[count],
                        "all_pool_calls": int(call.sum()),
                        "scored_calls": int(call[labelled].sum()),
                        "scored_positives_called": int(
                            (call[labelled] & (raw_y[labelled] == 1)).sum()
                        ),
                        "projected_mcc": observed,
                        "delta_vs_reference": observed - baseline,
                        "delta_ci025": float(np.quantile(deltas, 0.025)),
                        "delta_ci975": float(np.quantile(deltas, 0.975)),
                        "bootstrap_probability_positive": float((np.asarray(deltas) > 0).mean()),
                    }
                )
                bootstrap_rows.append(
                    {
                        "split_id": split,
                        "prior": prior,
                        "count_per_750": count,
                        "bootstrap_valid_replicates": len(deltas),
                    }
                )
            for held in np.unique(folds):
                training = labelled & (folds != held)
                ranked = [
                    (
                        projected_mcc(raw_y[training].astype(int), calls[c][training], prior),
                        c,
                    )
                    for c in counts
                ]
                inner_mcc, selected = max(ranked, key=lambda item: (item[0], -item[1]))
                nested_call[folds == held] = calls[selected][folds == held]
                held_labelled = labelled & (folds == held)
                choice_rows.append(
                    {
                        "split_id": split,
                        "prior": prior,
                        "held_fold": int(held),
                        "selected_count_per_750": selected,
                        "inner_projected_mcc": inner_mcc,
                        "held_scored_rows": int(held_labelled.sum()),
                        "held_selected_calls_scored": int(calls[selected][held_labelled].sum()),
                        "held_reference_calls_scored": int(reference_calls[held_labelled].sum()),
                        "held_selected_mcc": projected_mcc(
                            raw_y[held_labelled].astype(int), calls[selected][held_labelled], prior
                        ),
                        "held_reference_mcc": projected_mcc(
                            raw_y[held_labelled].astype(int), reference_calls[held_labelled], prior
                        ),
                    }
                )
            nested_mcc = projected_mcc(raw_y[labelled].astype(int), nested_call[labelled], prior)
            reference_mcc = projected_mcc(
                raw_y[labelled].astype(int), reference_calls[labelled], prior
            )
            nested_deltas = []
            for rows in samples:
                y = raw_y[rows].astype(int)
                delta = projected_mcc(y, nested_call[rows], prior) - projected_mcc(
                    y, reference_calls[rows], prior
                )
                if np.isfinite(delta):
                    nested_deltas.append(delta)
            policy_rows.append(
                {
                    "split_id": split,
                    "prior": prior,
                    "nested_projected_mcc": nested_mcc,
                    "reference_projected_mcc": reference_mcc,
                    "nested_delta_vs_reference": nested_mcc - reference_mcc,
                    "delta_ci025": float(np.quantile(nested_deltas, 0.025)),
                    "delta_ci975": float(np.quantile(nested_deltas, 0.975)),
                    "nested_calls_all": int(nested_call.sum()),
                    "nested_calls_scored": int(nested_call[labelled].sum()),
                }
            )
    curve = pl.DataFrame(curve_rows)
    choices = pl.DataFrame(choice_rows)
    curve.write_csv(outputs / "rate-curve.csv")
    choices.write_csv(outputs / "outer-fold-choices.csv")
    pl.DataFrame(policy_rows).write_csv(outputs / "nested-policy.csv")
    pl.DataFrame(bootstrap_rows).write_csv(outputs / "bootstrap-coverage.csv")
    summary = {
        "reference_count": reference,
        "counts_per_750": counts,
        "scored_population_priors": cfg["scored_population_priors"],
        "pool_rows_per_split": cfg["expected_pool_rows"],
        "labelled_rows_per_split": {
            split: load_oof(
                inputs / alias, cfg["model"], cfg["endpoint"], split, cfg["expected_pool_rows"]
            )["y_true"]
            .is_not_null()
            .sum()
            for split, alias in cfg["oof"].items()
        },
        "interpretation": cfg["protocol"],
    }
    (outputs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(curve.filter(pl.col("prior") == cfg["scored_population_priors"][1]), flush=True)


if __name__ == "__main__":
    main()
