import math

import pytest

from cyp_challenge import descriptors
from cyp_challenge.preprocessing import standardize


def test_salt_parent_deduplication():
    a = standardize("CC[NH+](CC)CCc1ccccc1.[Cl-]")
    b = standardize("CCN(CC)CCc1ccccc1")
    assert a["compound_id"] == b["compound_id"]
    assert a["status"] == "changed"


def test_stereoisomers_remain_distinct():
    a = standardize("C[C@H](O)C(=O)O")
    b = standardize("C[C@@H](O)C(=O)O")
    assert a["compound_id"] and b["compound_id"]
    assert a["compound_id"] != b["compound_id"]


@pytest.mark.parametrize("smiles", ["not-a-smiles", ""])
def test_invalid_input_is_quarantined(smiles):
    result = standardize(smiles)
    assert result["compound_id"] is None
    assert result["status"] == "invalid_input"


def test_standardization_idempotence():
    a = standardize("CC[NH+](CC)CCc1ccccc1.[Cl-]")
    b = standardize(a["standardized_smiles"])
    assert a["compound_id"] == b["compound_id"]
    assert b["status"] == "unchanged"


def test_full_descriptor_schema_and_explicit_missingness():
    descriptors.initialize()
    result = descriptors.calculate(("ethanol", "CCO"))
    assert set(result["rdkit"]) == {name for name, _ in descriptors.RDKIT}
    assert set(result["mordred"]) == {str(d) for d in descriptors.CALCULATOR.descriptors}
    for family in ("rdkit", "mordred"):
        for name, value in result[family].items():
            if value is None:
                assert f"{family}:{name}" in result["errors"]
            else:
                assert math.isfinite(value)
    assert result["rdkit"]["MolWt"] == pytest.approx(46.069, abs=0.001)


def test_excluded_metal_keeps_raw_identity():
    result = standardize("N.N.[Cl][Pt][Cl]")
    assert result["excluded"]
    assert result["compound_id"] is None
    assert result["raw_smiles"] == "N.N.[Cl][Pt][Cl]"


def test_exact_measurement_dedup_preserves_conflicts_and_test_ids(tmp_path, monkeypatch):
    import json
    from pathlib import Path

    import polars as pl

    from cyp_challenge.preprocessing import prepare

    inputs, output = tmp_path / "inputs", tmp_path / "outputs"
    inputs.mkdir()
    output.mkdir()
    monkeypatch.setenv("CYP_INPUT_DIR", str(inputs))
    monkeypatch.setenv("CYP_OUTPUT_DIR", str(output))
    monkeypatch.setenv("CYP_CONFIG_PATH", str(Path("configs/preprocessing/default.json").resolve()))
    endpoint = "CYP2D6_pIC50_direct_inhibition"
    specs = [
        ("TEST-BLINDED", [{"Molecule_Name": "test", "SMILES": "CCO"}]),
        ("TRAIN_Emax", [{"Molecule_Name": "train", "SMILES": "CCO", endpoint: 6.0}]),
        ("TRAIN_TDI", [{"Molecule_Name": "train", "SMILES": "CCO", endpoint: 5.0}]),
        ("TRAIN_inhibition", [{"Molecule_Name": "train", "SMILES": "CCO", endpoint: 5.0}] * 2),
        (
            "single-concentration-TRAIN",
            [
                {
                    "Molecule_Name": "aux",
                    "SMILES": "OCC",
                    "enzyme": "CYP2D6",
                    "log2fc_estimate": -1.0,
                }
            ],
        ),
    ]
    for index, (name, rows) in enumerate(specs):
        pl.DataFrame(rows).write_csv(inputs / f"{index:03d}-cyp-challenge-{name}.csv")
    prepare()
    summary = json.loads((output / "summary.json").read_text())
    assert summary["unique_compounds"] == 1
    assert summary["deduplicated_measurements"] == 2
    assert summary["conflicting_measurement_groups"] == 1
    assert summary["test_overlap_compounds"] == 1
    assert summary["files"]["cyp-challenge-TRAIN_inhibition.csv"]["exact_duplicates_removed"] == 1
    assert pl.read_parquet(output / "test_mapping.parquet")["molecule_name"].to_list() == ["test"]
    assert pl.read_parquet(output / "record_origins.parquet").height == 6
