from cyp_challenge.classification_preflight import (
    evaluate_preflight,
    exchangeable_effective_count,
)

THRESHOLDS = {
    "max_swap_fraction": 0.25,
    "max_mean_member_score_spearman": 0.95,
    "min_effective_vote_count": 2.0,
    "min_oof_swap_positive_rate_advantage": 0.02,
    "max_oof_swap_fisher_p_two_sided": 0.1,
    "max_test_minus_oof_member_spearman": 0.05,
    "label_neighborhood_veto_p": 0.05,
}


def test_exchangeable_effective_count():
    assert exchangeable_effective_count(6, 0.0) == 6
    assert abs(exchangeable_effective_count(6, 0.8646248799) - 1.127157577) < 1e-8


def test_support6_like_candidate_fails_all_safety_gates():
    metrics = {
        "swap_fraction": 27 / 49,
        "mean_test_member_score_spearman": 0.9700,
        "effective_vote_count": 1.1272,
        "oof_candidate_minus_control_swap_positive_rate": 0.02817,
        "oof_swap_fisher_p_two_sided": 0.8652,
        "test_minus_oof_member_spearman": 0.1138,
        "control_minus_candidate_top10_weighted_positive_rate": 0.1412,
        "top10_weighted_positive_rate_permutation_p": 0.0003,
    }
    status, checks = evaluate_preflight(metrics, THRESHOLDS)
    assert status == "FAIL"
    assert {row["check"] for row in checks if row["status"] == "FAIL"} == {
        "swap_fraction",
        "member_score_redundancy",
        "effective_vote_count",
        "oof_swap_evidence",
        "deployment_correlation_inflation",
        "label_neighborhood_veto",
    }


def test_well_supported_candidate_passes():
    metrics = {
        "swap_fraction": 0.12,
        "mean_test_member_score_spearman": 0.82,
        "effective_vote_count": 2.8,
        "oof_candidate_minus_control_swap_positive_rate": 0.12,
        "oof_swap_fisher_p_two_sided": 0.03,
        "test_minus_oof_member_spearman": 0.02,
        "control_minus_candidate_top10_weighted_positive_rate": -0.04,
        "top10_weighted_positive_rate_permutation_p": 0.02,
    }
    status, checks = evaluate_preflight(metrics, THRESHOLDS)
    assert status == "PASS"
    assert all(row["status"] == "PASS" for row in checks)
