"""Compose the submitted A4 equal-seven score on the complete UMAP pool."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def main() -> None:
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    selection = pl.read_parquet(inputs / cfg["selection"]).filter(pl.col("model") == cfg["model"])
    if selection.height != 1:
        raise ValueError("Exact equal-seven selection missing")
    choices = json.loads(selection["members"][0])
    full_choice = next(row["members"] for row in choices if row["held_fold"] == -1)
    for fold in range(5):
        if next(row["members"] for row in choices if row["held_fold"] == fold) != full_choice:
            raise ValueError("Outer-fold member selection differs from submitted recipe")
    if len(full_choice) != 7 or len(set(full_choice)) != 7:
        raise ValueError("Submitted member selection is not seven unique members")
    bank = pl.concat(
        [pl.read_parquet(inputs / cfg[key]) for key in ("historical", "jeremy", "trees")]
    )
    bank = bank.filter(
        (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("split_id") == cfg["split_id"])
    )
    if sorted(bank["model"].unique().to_list()) != sorted(full_choice):
        raise ValueError("Member bank differs from exact submitted selection")
    frame = None
    for member in full_choice:
        part = bank.filter(pl.col("model") == member).select(
            "compound_id", "fold", "group", "y_true", pl.col("score").alias(member)
        )
        if (
            part.height != cfg["expected_pool_rows"]
            or part["compound_id"].n_unique() != part.height
        ):
            raise ValueError(f"Incomplete member pool: {member}")
        frame = (
            part
            if frame is None
            else frame.join(
                part.select("compound_id", pl.col(member)), on="compound_id", validate="1:1"
            )
        )
    frame = frame.sort("compound_id")
    folds = frame["fold"].to_numpy()
    scored = frame["y_true"].is_not_null().to_numpy()
    rank_all = []
    rank_labelled = []
    for member in full_choice:
        values = frame[member].to_numpy()
        full_rank = np.empty(frame.height, dtype=float)
        labelled_rank = np.full(frame.height, np.nan)
        for fold in range(5):
            mask = folds == fold
            full_rank[mask] = unit_rank(values[mask])
            labelled_mask = mask & scored
            labelled_rank[labelled_mask] = unit_rank(values[labelled_mask])
        rank_all.append(full_rank)
        rank_labelled.append(labelled_rank)
    full_score = np.mean(rank_all, axis=0)
    labelled_score = np.nanmean(np.column_stack(rank_labelled)[scored], axis=1)
    original = pl.read_parquet(inputs / cfg["saved_labelled_oof"]).filter(
        pl.col("model") == cfg["model"]
    )
    compared = (
        frame.filter(pl.Series(scored))
        .select("compound_id")
        .with_columns(pl.Series("replayed", labelled_score))
        .join(
            original.select("compound_id", pl.col("prediction").alias("original")),
            on="compound_id",
            validate="1:1",
        )
    )
    if compared.height != int(scored.sum()):
        raise ValueError("Labelled composition audit lacks rows")
    labelled_difference = np.abs(compared["replayed"].to_numpy() - compared["original"].to_numpy())
    test = pl.read_parquet(inputs / cfg["saved_test"]).filter(pl.col("model") == cfg["model"])
    submitted = pl.read_parquet(inputs / cfg["submitted_test"]).select(
        "compound_id", "a4_prediction", "a4_call"
    )
    test_audit = test.join(submitted, on="compound_id", validate="1:1")
    if test_audit.height != cfg["expected_test_rows"]:
        raise ValueError("Blind test score audit lacks rows")
    if not test_audit["call"].equals(test_audit["a4_call"]):
        raise ValueError("Saved equal-seven calls differ from accepted submission")
    if not np.allclose(test_audit["prediction"], test_audit["a4_prediction"], atol=0, rtol=0):
        raise ValueError("Saved equal-seven scores differ from accepted package")
    if test_audit["call"].sum() != cfg["reference_count"]:
        raise ValueError("Wrong submitted A4 positive count")
    prediction = frame.select("compound_id", "fold", "group", "y_true").with_columns(
        pl.lit(cfg["endpoint"]).alias("endpoint"),
        pl.lit(cfg["split_id"]).alias("split_id"),
        pl.lit("submitted_equal7_full_pool").alias("model"),
        pl.Series("score", full_score),
    )
    prediction.write_parquet(outputs / "all-pool-oof.parquet")
    audit = {
        "members": full_choice,
        "pool_rows": frame.height,
        "labelled_rows": int(scored.sum()),
        "labelled_recomposition_max_abs_diff": float(labelled_difference.max()),
        "labelled_recomposition_mae": float(labelled_difference.mean()),
        "labelled_recomposition_spearman": float(
            spearmanr(compared["replayed"], compared["original"]).statistic
        ),
        "submitted_test_calls_exact": True,
        "submitted_test_count": int(test_audit["call"].sum()),
    }
    (outputs / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
