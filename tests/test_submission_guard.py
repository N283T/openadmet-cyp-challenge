"""Guard deliberate payload resubmission without weakening uncertain-delivery safety."""

import pytest

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
