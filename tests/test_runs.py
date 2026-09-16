from pathlib import Path

import pytest

from cyp_challenge.runs import resolve_inputs


def test_resolve_inputs_prefers_explicit_aliases(tmp_path):
    source = tmp_path / "000-oof.parquet"
    source.write_text("value")

    resolved = resolve_inputs([f"oof.parquet={source}"])

    assert resolved == [(source.resolve(), "oof.parquet")]


def test_resolve_inputs_keeps_legacy_numeric_prefix(tmp_path):
    source = tmp_path / "oof.parquet"
    source.write_text("value")

    resolved = resolve_inputs([Path(source)])

    assert resolved == [(source.resolve(), "000-oof.parquet")]


def test_resolve_inputs_rejects_invalid_or_duplicate_aliases(tmp_path):
    first = tmp_path / "first.parquet"
    second = tmp_path / "second.parquet"
    first.write_text("first")
    second.write_text("second")

    with pytest.raises(ValueError, match="Invalid input alias"):
        resolve_inputs([f"nested/oof.parquet={first}"])
    with pytest.raises(ValueError, match="Duplicate input snapshot name"):
        resolve_inputs([f"oof.parquet={first}", f"oof.parquet={second}"])


def test_resolve_inputs_checks_paths_before_run_allocation(tmp_path):
    missing = tmp_path / "missing.parquet"

    with pytest.raises(FileNotFoundError):
        resolve_inputs([f"oof.parquet={missing}"])
