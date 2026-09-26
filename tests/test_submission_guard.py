"""Guard deliberate payload resubmission without weakening uncertain-delivery safety."""

from datetime import UTC, datetime, timedelta

import polars as pl
import pytest

from cyp_challenge import submission
from cyp_challenge.submission import check_previous_attempts


def test_sent_payload_requires_new_id_and_explicit_flag():
    previous = [("old", "sent")]
    with pytest.raises(ValueError, match="--resubmit-sent"):
        check_previous_attempts(previous, "new", False, False)
    check_previous_attempts(previous, "new", True, False)
    with pytest.raises(ValueError, match="already been attempted"):
        check_previous_attempts(previous, "old", True, False)


def test_unknown_delivery_blocks_even_explicit_resubmission():
    with pytest.raises(ValueError, match="uncertain"):
        check_previous_attempts([("old", "unknown")], "new", True, True)


def test_rejected_payload_still_requires_its_existing_retry_flag():
    with pytest.raises(ValueError, match="--retry-rejected"):
        check_previous_attempts([("old", "rejected")], "new", False, False)
    check_previous_attempts([("old", "rejected")], "new", False, True)


def test_target_changes_reports_each_endpoint(tmp_path):
    endpoints = submission.ENDPOINTS["classification"]
    old = pl.DataFrame(
        {
            "Molecule_Name": ["a", "b"],
            "SMILES": ["C", "CC"],
            endpoints[0]: [True, False],
            endpoints[1]: [False, True],
        }
    )
    new = old.with_columns(pl.Series(endpoints[0], [False, False]))
    old_path, new_path = tmp_path / "old.parquet", tmp_path / "new.parquet"
    old.write_parquet(old_path)
    new.write_parquet(new_path)
    assert submission.target_changes(new_path, old_path, "classification") == {
        endpoints[0]: 1,
        endpoints[1]: 0,
    }


def test_live_confirmation_requires_tty_and_exact_id(monkeypatch):
    monkeypatch.setattr(submission.sys.stdin, "isatty", lambda: False)
    with pytest.raises(ValueError, match="interactive terminal"):
        submission.confirm_live_submission("prepared-id", "regression", None, None)
    monkeypatch.setattr(submission.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: "SUBMIT another-id")
    with pytest.raises(ValueError, match="did not match"):
        submission.confirm_live_submission("prepared-id", "regression", None, None)
    monkeypatch.setattr("builtins.input", lambda _: "SUBMIT prepared-id")
    submission.confirm_live_submission("prepared-id", "regression", None, None)


def test_cooldown_blocks_until_twelve_hours_have_elapsed():
    sent_at = datetime(2026, 9, 27, tzinfo=UTC)
    with pytest.raises(ValueError, match="12-hour track cooldown"):
        submission.check_cooldown(sent_at, sent_at + timedelta(hours=12, seconds=-1))
    assert submission.check_cooldown(sent_at, sent_at + timedelta(hours=12)) == (
        sent_at + timedelta(hours=12)
    )
