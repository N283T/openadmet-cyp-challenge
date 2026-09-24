import numpy as np
import polars as pl
import pytest

from cyp_challenge.classification_evaluation import (
    evaluate_rates,
    fixed_rate_calls,
    prepare_predictions,
)


def test_full_pool_selection_counts_unscored_compounds_in_call_budget():
    rows = []
    assignments = []
    for fold in range(3):
        for index, (label, control, candidate) in enumerate(
            [(1, 0.7, 0.8), (0, 0.6, 0.7), (None, 0.9, 0.9), (None, 0.1, 0.1)]
        ):
            compound = f"{fold}-{index}"
            assignments.append(
                {"strategy": "umap", "compound_id": compound, "fold": fold, "group": compound}
            )
            for model, score in [("control", control), ("candidate", candidate)]:
                rows.append(
                    {"compound_id": compound, "model": model, "score": score, "y_true": label}
                )
    data = prepare_predictions(
        pl.DataFrame(rows), pl.DataFrame(assignments), ["control", "candidate"], "umap"
    )
    summary, _, decisions, calls = evaluate_rates(
        data, ["control", "candidate"], [0.25], 0.1, 12, 100, 42
    )
    assert summary["deployment_decision_supported"]
    assert summary["unscored_rows"] == 6
    assert decisions.height == 9
    assert calls.filter(pl.col("compound_id").str.ends_with("-2"))["control_call"].all()
    assert not calls.filter(pl.col("compound_id").str.ends_with("-0"))["control_call"].any()


def test_rate_selection_uses_only_other_folds():
    # The third fold favors a wider rate, but it must not choose its own rate.
    rows = []
    assignments = []
    labels = [[1, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0]]
    for fold in range(3):
        for index in range(4):
            compound = f"{fold}-{index}"
            assignments.append(
                {"strategy": "umap", "compound_id": compound, "fold": fold, "group": compound}
            )
            for model in ["control", "candidate"]:
                rows.append(
                    {
                        "compound_id": compound,
                        "model": model,
                        "score": 4 - index,
                        "y_true": labels[fold][index],
                    }
                )
    data = prepare_predictions(
        pl.DataFrame(rows), pl.DataFrame(assignments), ["control", "candidate"], "umap"
    )
    _, _, decisions, _ = evaluate_rates(
        data, ["control", "candidate"], [0.25, 0.5], 0.25, 12, 30, 7
    )
    selected = decisions.filter((pl.col("model") == "control") & (pl.col("held_fold") == 2))
    assert selected["selected_rate"].item() == 0.25


def test_rejects_partial_label_disagreement_and_bad_rate():
    predictions = pl.DataFrame(
        {
            "compound_id": ["a", "a"],
            "model": ["control", "candidate"],
            "score": [0.1, 0.2],
            "y_true": [1, None],
        }
    )
    assignments = pl.DataFrame(
        {"strategy": ["umap"], "compound_id": ["a"], "fold": [0], "group": ["a"]}
    )
    with pytest.raises(ValueError, match="different molecule coverage or labels"):
        prepare_predictions(predictions, assignments, ["control", "candidate"], "umap")
    with pytest.raises(ValueError, match="rate"):
        fixed_rate_calls(np.array([0.1, 0.2]), np.array([0, 0]), 1.0)


def test_incomplete_pool_cannot_support_deployment_decision():
    rows = []
    assignments = []
    for fold in range(3):
        for index, y in enumerate([0, 1]):
            compound_id = f"{fold}-{index}"
            assignments.append(
                {"strategy": "umap", "compound_id": compound_id, "fold": fold, "group": compound_id}
            )
            for model in ["control", "candidate"]:
                rows.append(
                    {"compound_id": compound_id, "model": model, "score": float(index), "y_true": y}
                )
    data = prepare_predictions(
        pl.DataFrame(rows), pl.DataFrame(assignments), ["control", "candidate"], "umap"
    )
    with pytest.raises(ValueError, match="Incomplete molecular pool"):
        evaluate_rates(data, ["control", "candidate"], [0.5], 0.1, 9, 10, 42)
    summary, *_ = evaluate_rates(
        data,
        ["control", "candidate"],
        [0.5],
        0.1,
        9,
        10,
        42,
        allow_partial_pool=True,
    )
    assert not summary["deployment_decision_supported"]
