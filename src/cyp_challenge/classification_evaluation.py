"""Decision-level classification evaluation with a selectable positive-call rate."""

from __future__ import annotations

import numpy as np
import polars as pl


def prepare_predictions(
    predictions: pl.DataFrame,
    assignments: pl.DataFrame,
    models: list[str],
    strategy: str,
) -> pl.DataFrame:
    """Align all model scores to one molecular pool and its saved outer folds."""
    required = {"compound_id", "model", "score", "y_true"}
    if not required <= set(predictions.columns):
        raise ValueError(f"Predictions missing {sorted(required - set(predictions.columns))}")
    if len(models) != 2 or len(set(models)) != 2:
        raise ValueError("Exactly two distinct control/candidate models are required")
    selected = predictions.filter(pl.col("model").is_in(models)).select(
        "compound_id", "model", "score", "y_true"
    )
    if selected.is_empty() or selected.select("compound_id", "model").is_duplicated().any():
        raise ValueError("Predictions must have unique compound/model keys")
    if selected["score"].null_count() or not np.isfinite(selected["score"].to_numpy()).all():
        raise ValueError("Every prediction score must be finite")
    membership = assignments.filter(pl.col("strategy") == strategy).select(
        "compound_id", "fold", "group"
    )
    if membership.is_empty() or membership.select("compound_id").is_duplicated().any():
        raise ValueError("Missing or duplicated split assignments")
    joined = selected.join(membership, on="compound_id", how="left", validate="m:1")
    if joined["fold"].null_count() or joined["group"].null_count():
        raise ValueError("A predicted molecule has no split assignment")
    counts = joined.group_by("compound_id").agg(
        pl.col("model").n_unique().alias("n_models"),
        pl.col("y_true").drop_nulls().n_unique().alias("n_labels"),
        pl.col("y_true").is_not_null().sum().alias("n_scored"),
    )
    if (
        (counts["n_models"] != len(models)).any()
        or (counts["n_labels"] > 1).any()
        or (~counts["n_scored"].is_in([0, len(models)])).any()
    ):
        raise ValueError("Models have different molecule coverage or labels")
    if not set(joined["y_true"].drop_nulls().unique().to_list()) <= {0, 1}:
        raise ValueError("Labels must be binary or null for unscored molecules")
    if joined["fold"].n_unique() < 3:
        raise ValueError("At least three outer folds are required")
    return joined.sort("model", "compound_id")


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    """Choose the highest scores from every full outer-fold pool."""
    if not 0 < rate < 1:
        raise ValueError("Positive-call rate must be between zero and one")
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, min(len(rows) - 1, int(round(rate * len(rows)))))
        calls[rows[np.argsort(-scores[rows], kind="stable")[:count]]] = True
    return calls


def projected_mcc(y: np.ndarray, calls: np.ndarray, prior: float) -> float:
    """Apply a prespecified scored-population prior to observed TPR and FPR."""
    if not 0 < prior < 1 or not (y == 1).any() or not (y == 0).any():
        return float("nan")
    tpr = float(calls[y == 1].mean())
    fpr = float(calls[y == 0].mean())
    tp, fn = prior * tpr, prior * (1 - tpr)
    fp, tn = (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def evaluate_rates(
    data: pl.DataFrame,
    models: list[str],
    rates: list[float],
    prior: float,
    expected_pool_size: int,
    bootstrap_replicates: int,
    bootstrap_seed: int,
    *,
    allow_partial_pool: bool = False,
) -> tuple[dict, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Select each model's call rate without the reported outer fold."""
    if not rates or len(set(rates)) != len(rates) or any(not 0 < r < 1 for r in rates):
        raise ValueError("Rates must be distinct values strictly between zero and one")
    reference = data.filter(pl.col("model") == models[0]).sort("compound_id")
    ids = reference["compound_id"].to_list()
    folds = reference["fold"].to_numpy()
    groups = reference["group"].to_numpy()
    y_all = reference["y_true"].to_numpy()
    scored = ~np.isnan(y_all.astype(float))
    y = y_all[scored].astype(int)
    if len(ids) != expected_pool_size and not allow_partial_pool:
        raise ValueError(
            f"Incomplete molecular pool: {len(ids)} predictions for "
            f"{expected_pool_size} assigned molecules"
        )
    if expected_pool_size < len(ids):
        raise ValueError("Expected molecular pool cannot be smaller than predictions")
    curves: list[dict] = []
    decisions: list[dict] = []
    selected_calls: dict[str, np.ndarray] = {}
    selected_rates: dict[str, list[float]] = {}
    call_cache: dict[str, dict[float, np.ndarray]] = {}
    all_folds = np.unique(folds)
    for model in models:
        part = data.filter(pl.col("model") == model).sort("compound_id")
        if part["compound_id"].to_list() != ids:
            raise ValueError("Model alignment failed")
        scores = part["score"].to_numpy()
        model_calls: dict[float, np.ndarray] = {}
        for rate in rates:
            calls = fixed_rate_calls(scores, folds, rate)
            model_calls[rate] = calls
            for fold in all_folds:
                mask = scored & (folds == fold)
                curves.append(
                    {
                        "model": model,
                        "rate": rate,
                        "fold": int(fold),
                        "pool_rows": int((folds == fold).sum()),
                        "scored_rows": int(mask.sum()),
                        "positive_calls_all": int(calls[folds == fold].sum()),
                        "positive_calls_scored": int(calls[mask].sum()),
                        "projected_mcc": projected_mcc(y_all[mask].astype(int), calls[mask], prior),
                    }
                )
        call_cache[model] = model_calls
        nested = np.zeros(len(ids), dtype=bool)
        choices = []
        for held in all_folds:
            train = scored & (folds != held)
            ranking = [
                (projected_mcc(y_all[train].astype(int), model_calls[r][train], prior), r)
                for r in rates
            ]
            best_value, best_rate = max(ranking, key=lambda pair: (pair[0], -pair[1]))
            valid = folds == held
            nested[valid] = model_calls[best_rate][valid]
            choices.append(best_rate)
            decisions.append(
                {
                    "model": model,
                    "held_fold": int(held),
                    "selected_rate": best_rate,
                    "inner_projected_mcc": best_value,
                    "outer_scored_calls": int((nested & valid & scored).sum()),
                    "outer_scored_positives": int((nested & valid & scored & (y_all == 1)).sum()),
                    "selected_model": model,
                }
            )
        selected_calls[model] = nested
        selected_rates[model] = choices
    selector_calls = np.zeros(len(ids), dtype=bool)
    selector_choices = []
    for held in all_folds:
        train = scored & (folds != held)
        options = [
            (
                projected_mcc(y_all[train].astype(int), call_cache[model][rate][train], prior),
                model,
                rate,
            )
            for model in models
            for rate in rates
        ]
        best_value, best_model, best_rate = max(
            options, key=lambda option: (option[0], -models.index(option[1]), -option[2])
        )
        valid = folds == held
        selector_calls[valid] = call_cache[best_model][best_rate][valid]
        selector_choices.append({"fold": int(held), "model": best_model, "rate": best_rate})
        decisions.append(
            {
                "model": "selected_policy",
                "held_fold": int(held),
                "selected_rate": best_rate,
                "inner_projected_mcc": best_value,
                "outer_scored_calls": int((selector_calls & valid & scored).sum()),
                "outer_scored_positives": int(
                    (selector_calls & valid & scored & (y_all == 1)).sum()
                ),
                "selected_model": best_model,
            }
        )
    selected_calls["selected_policy"] = selector_calls
    evaluation = []
    for model in [*models, "selected_policy"]:
        calls = selected_calls[model]
        evaluation.append(
            {
                "model": model,
                "projected_mcc": projected_mcc(y, calls[scored], prior),
                "total_calls": int(calls.sum()),
                "scored_calls": int(calls[scored].sum()),
                "scored_true_positives": int((calls[scored] & (y == 1)).sum()),
                "selected_rates": selected_rates.get(model),
            }
        )
    unique_groups = np.unique(groups[scored])
    rows_by_group = {group: np.flatnonzero(scored & (groups == group)) for group in unique_groups}
    rng = np.random.default_rng(bootstrap_seed)
    deltas = {model: [] for model in [models[1], "selected_policy"]}
    for _ in range(bootstrap_replicates):
        sample = rng.choice(unique_groups, len(unique_groups), replace=True)
        rows = np.concatenate([rows_by_group[group] for group in sample])
        truth = y_all[rows].astype(int)
        control_score = projected_mcc(truth, selected_calls[models[0]][rows], prior)
        for model in deltas:
            difference = projected_mcc(truth, selected_calls[model][rows], prior) - control_score
            if np.isfinite(difference):
                deltas[model].append(difference)
    if any(not value for value in deltas.values()):
        raise ValueError("Bootstrap produced no valid paired comparisons")

    def delta_summary(model: str, evaluation_index: int) -> dict:
        delta = np.asarray(deltas[model])
        return {
            "point": evaluation[evaluation_index]["projected_mcc"] - evaluation[0]["projected_mcc"],
            "bootstrap_mean": float(delta.mean()),
            "ci025": float(np.quantile(delta, 0.025)),
            "ci975": float(np.quantile(delta, 0.975)),
            "probability_positive": float((delta > 0).mean()),
            "valid_replicates": len(delta),
        }

    summary = {
        "pool_rows": len(ids),
        "expected_pool_rows": expected_pool_size,
        "pool_coverage": len(ids) / expected_pool_size,
        "scored_rows": int(scored.sum()),
        "unscored_rows": int((~scored).sum()),
        "target_prior": prior,
        "rates": rates,
        "models": models,
        "nested_evaluation": evaluation,
        "selector_choices": selector_choices,
        "paired_delta_candidate_minus_control": delta_summary(models[1], 1),
        "paired_delta_selector_minus_control": delta_summary("selected_policy", 2),
        "deployment_decision_supported": len(ids) == expected_pool_size,
    }
    return (
        summary,
        pl.DataFrame(curves),
        pl.DataFrame(decisions),
        pl.DataFrame(
            {
                "compound_id": ids,
                "fold": folds,
                "y_true": y_all,
                **{f"{model}_call": selected_calls[model] for model in models},
            }
        ),
    )
