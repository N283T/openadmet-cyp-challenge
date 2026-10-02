"""Exact feasible MCC-difference bounds with missing outcomes, never imputed truth."""

import json
import os
from itertools import product
from pathlib import Path

import numpy as np
import polars as pl
from review_cyp_d6_causal_probe import calls_at_rate
from sklearn.metrics import matthews_corrcoef


def bounds(y, known, before, after):
    n = len(y)
    ka, kb = int(before.sum()), int(after.sum())
    fa = 1 / np.sqrt(ka * (n - ka)) if 0 < ka < n else 0.0
    fb = 1 / np.sqrt(kb * (n - kb)) if 0 < kb < n else 0.0
    fixed_positive = int(y[known].sum())
    fixed_contrast = fb * np.sum(y[known] * after[known]) - fa * np.sum(y[known] * before[known])
    coefficients = np.sort((fb * after - fa * before)[~known])
    positive = fixed_positive + np.arange(len(coefficients) + 1)
    denominator = np.sqrt(positive * (n - positive))
    sums_low = np.r_[0.0, np.cumsum(coefficients)]
    sums_high = np.r_[0.0, np.cumsum(coefficients[::-1])]
    low_num = n * (fixed_contrast + sums_low) - positive * (kb * fb - ka * fa)
    high_num = n * (fixed_contrast + sums_high) - positive * (kb * fb - ka * fa)
    low = np.divide(low_num, denominator, out=np.zeros_like(low_num), where=denominator > 0)
    high = np.divide(high_num, denominator, out=np.zeros_like(high_num), where=denominator > 0)
    lower_idx, upper_idx = int(np.argmin(low)), int(np.argmax(high))
    known_tp_delta = int(np.sum(y[known] * (after[known].astype(int) - before[known].astype(int))))
    removed_unknown = int(((before & ~after) & ~known).sum())
    added_unknown = int(((after & ~before) & ~known).sum())
    return dict(
        population_n=n,
        known_n=int(known.sum()),
        unknown_n=int((~known).sum()),
        before_calls=ka,
        after_calls=kb,
        known_tp_delta=known_tp_delta,
        removed_unknown=removed_unknown,
        added_unknown=added_unknown,
        lower_delta=float(low[lower_idx]),
        upper_delta=float(high[upper_idx]),
        lower_unknown_positives=lower_idx,
        upper_unknown_positives=upper_idx,
        possible_negative=bool(low.min() < 0),
        possible_positive=bool(high.max() > 0),
        equal_count_reversal_witness=(
            known_tp_delta + 1
            if ka == kb and known_tp_delta >= 0 and removed_unknown > known_tp_delta
            else None
        ),
    )


def self_check():
    """Compare the analytic feasible interval to exhaustive, independent tiny enumerations."""
    rng = np.random.default_rng(20261002)
    for number in range(24):
        y = rng.integers(0, 2, 8).astype(float)
        known = np.arange(8) < 3
        before, after = rng.integers(0, 2, (2, 8)).astype(bool)
        if number == 0:
            before[:] = False
        if number == 1:
            after[:] = True
        result = bounds(y, known, before, after)
        exact = []
        for assignment in product([0, 1], repeat=5):
            truth = y.copy()
            truth[~known] = assignment
            exact.append(matthews_corrcoef(truth, after) - matthews_corrcoef(truth, before))
        np.testing.assert_allclose(result["lower_delta"], min(exact), atol=1e-12)
        np.testing.assert_allclose(result["upper_delta"], max(exact), atol=1e-12)
    return 24 * 32


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    exhaustive_checked = self_check()
    rows = []
    for case in cfg["cases"]:
        raw = pl.read_parquet(inp / case["before_file"])
        second = pl.read_parquet(inp / case["after_file"])
        if case["kind"] == "a4_control":
            source = pl.read_parquet(inp / case["source"]).select("compound_id", "selected")
            for split in ["umap", "scaffold"]:
                for seed in cfg["seeds"]:
                    parts = []
                    for data, model in [(raw, case["before_model"]), (second, case["after_model"])]:
                        part = (
                            data.filter(
                                (pl.col("split_id") == split)
                                & (pl.col("seed") == seed)
                                & (pl.col("model") == model)
                                & (pl.col("policy") == "centered_fold_mean")
                                & pl.col("label").is_finite()
                            )
                            .join(source, on="compound_id", validate="1:1")
                            .sort("compound_id")
                        )
                        parts.append(part)
                    a, b = parts
                    assert a["compound_id"].equals(b["compound_id"])
                    y = a["label"].to_numpy()
                    known = a["selected"].to_numpy()
                    result = bounds(y, known, a["call"].to_numpy(), b["call"].to_numpy())
                    actual = float(
                        matthews_corrcoef(y, b["call"]) - matthews_corrcoef(y, a["call"])
                    )
                    assert result["lower_delta"] - 1e-12 <= actual <= result["upper_delta"] + 1e-12
                    rows.append(
                        dict(
                            case=case["name"],
                            split_id=split,
                            seed=seed,
                            scope="known A4 population; weak classes hidden only for bounds",
                            transport_calibration_eligible=True,
                            actual_delta=actual,
                            **result,
                        )
                    )
        else:
            parts = []
            for data, model in [(raw, case["before_model"]), (second, case["after_model"])]:
                model_col = "method" if "method" in data.columns else "model"
                score_col = "prediction" if "prediction" in data.columns else "score"
                parts.append(
                    data.filter(pl.col(model_col) == model)
                    .select("compound_id", "fold", "y_true", pl.col(score_col).alias("score"))
                    .sort("compound_id")
                )
            a, b = parts
            assert a["compound_id"].equals(b["compound_id"]) and a["fold"].equals(b["fold"])
            np.testing.assert_allclose(a["y_true"], b["y_true"], equal_nan=True)
            y = a["y_true"].to_numpy()
            known = np.isfinite(y)
            ca = calls_at_rate(a["score"].to_numpy(), a["fold"].to_numpy(), case["before_rate"])
            cb = calls_at_rate(b["score"].to_numpy(), b["fold"].to_numpy(), case["after_rate"])
            rows.append(
                dict(
                    case=case["name"],
                    split_id=case["split"],
                    seed=None,
                    scope=case["scope"],
                    transport_calibration_eligible=case["complete_pool"],
                    actual_delta=None,
                    **bounds(y, known, ca, cb),
                )
            )
    pl.DataFrame(rows).write_csv(out / "partial-label-bounds.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                results=rows,
                no_unknown_truth_assigned=True,
                no_upload=True,
                exhaustive_checked_assignments=exhaustive_checked,
                limits="Exact bounds over all unknown binary outcomes "
                "on the stated cached prediction "
                "population. They bound that population, not the separate hidden Live batch. "
                "Labelled-only historical successes cannot calibrate "
                "a complete-population selector. "
                "Witness counts describe feasible counterexamples, "
                "never reconstructed hidden labels.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
