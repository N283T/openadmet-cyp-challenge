"""Reusable decision rules for classification candidate preflight checks."""

from __future__ import annotations

from typing import Any


def exchangeable_effective_count(member_count: int, mean_correlation: float) -> float:
    """Approximate the independent vote count under exchangeable correlation."""
    if member_count < 1:
        raise ValueError("member_count must be positive")
    denominator = 1 + (member_count - 1) * mean_correlation
    if denominator <= 0:
        raise ValueError("correlation produces a non-positive denominator")
    return member_count / denominator


def evaluate_preflight(
    metrics: dict[str, float | int], thresholds: dict[str, float]
) -> tuple[str, list[dict[str, Any]]]:
    """Evaluate fixed candidate-safety gates and return overall status and rows."""
    checks: list[dict[str, Any]] = []

    def add(name: str, passed: bool, value: float, threshold: str, detail: str) -> None:
        checks.append(
            {
                "check": name,
                "status": "PASS" if passed else "FAIL",
                "value": value,
                "threshold": threshold,
                "detail": detail,
            }
        )

    swap_fraction = float(metrics["swap_fraction"])
    max_swap = float(thresholds["max_swap_fraction"])
    add(
        "swap_fraction",
        swap_fraction <= max_swap,
        swap_fraction,
        f"<= {max_swap}",
        "Fraction of control positive calls replaced by the candidate.",
    )

    member_spearman = float(metrics["mean_test_member_score_spearman"])
    max_spearman = float(thresholds["max_mean_member_score_spearman"])
    add(
        "member_score_redundancy",
        member_spearman <= max_spearman,
        member_spearman,
        f"<= {max_spearman}",
        "Mean pairwise blind-test rank correlation among claimed ensemble members.",
    )

    effective_votes = float(metrics["effective_vote_count"])
    min_votes = float(thresholds["min_effective_vote_count"])
    add(
        "effective_vote_count",
        effective_votes >= min_votes,
        effective_votes,
        f">= {min_votes}",
        "Exchangeable-correlation approximation from member call phi.",
    )

    oof_advantage = float(metrics["oof_candidate_minus_control_swap_positive_rate"])
    min_advantage = float(thresholds["min_oof_swap_positive_rate_advantage"])
    oof_p = float(metrics["oof_swap_fisher_p_two_sided"])
    max_oof_p = float(thresholds["max_oof_swap_fisher_p_two_sided"])
    add(
        "oof_swap_evidence",
        oof_advantage >= min_advantage and oof_p <= max_oof_p,
        oof_p,
        f"advantage >= {min_advantage} and p <= {max_oof_p}",
        f"Observed candidate-minus-control swap positive-rate advantage={oof_advantage:.6f}.",
    )

    inflation = float(metrics["test_minus_oof_member_spearman"])
    max_inflation = float(thresholds["max_test_minus_oof_member_spearman"])
    add(
        "deployment_correlation_inflation",
        inflation <= max_inflation,
        inflation,
        f"<= {max_inflation}",
        "Increase in member score correlation from honest OOF to blind test.",
    )

    neighborhood_delta = float(
        metrics["control_minus_candidate_top10_weighted_positive_rate"]
    )
    neighborhood_p = float(metrics["top10_weighted_positive_rate_permutation_p"])
    veto_p = float(thresholds["label_neighborhood_veto_p"])
    neighborhood_pass = not (neighborhood_delta > 0 and neighborhood_p <= veto_p)
    add(
        "label_neighborhood_veto",
        neighborhood_pass,
        neighborhood_delta,
        f"not (control-candidate > 0 and p <= {veto_p})",
        f"Permutation p={neighborhood_p:.6f}; positive values favor removed control calls.",
    )

    return ("PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL", checks)
