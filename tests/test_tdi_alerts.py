import pytest

from cyp_challenge.tdi_alerts import TDI_ALERTS, alert_matrix, compiled_alerts


def test_all_tdi_alert_smarts_compile():
    assert len(compiled_alerts()) == len(TDI_ALERTS) == 12


@pytest.mark.parametrize(
    ("name", "smiles"),
    [
        ("terminal_alkyne", "CC#C"),
        ("terminal_alkene", "C=CC"),
        ("furan", "c1ccoc1"),
        ("thiophene", "c1ccsc1"),
        ("methylenedioxyphenyl", "c1ccc2c(c1)OCO2"),
        ("primary_secondary_amine", "CCNCC"),
        ("hydrazine", "CNN"),
        ("aniline_aminophenol", "Nc1ccccc1"),
        ("thioamide", "CC(=S)N"),
        ("thiourea", "NC(=S)N"),
        ("cyclopropylamine", "NC1CC1"),
        ("n_dealkylatable_tertiary_amine", "CN(C)CC"),
    ],
)
def test_expected_example_fires(name, smiles):
    matrix, names = alert_matrix([smiles])
    assert matrix[0, names.index(name)] >= 1


def test_aggregates_match_individual_alert_counts():
    matrix, names = alert_matrix(["CC", "CN(C)CC"])
    alert_count = len(TDI_ALERTS)
    assert matrix[0, names.index("any_tdi_alert")] == 0
    assert matrix[1, names.index("any_tdi_alert")] == 1
    assert matrix[1, names.index("tdi_alert_count")] == matrix[1, :alert_count].sum()
