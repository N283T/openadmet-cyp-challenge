import json

import pytest

from cyp_challenge.submission import load_profile, response_status


def test_response_status_does_not_infer_delivery_from_ambiguous_text():
    assert (
        response_status(
            {"value": "Submission received from 'user' for the Regression Prediction track."}
        )
        == "sent"
    )
    assert response_status({"value": "Error: Please wait before submitting again."}) == "rejected"
    assert response_status("Please wait while processing") == "unknown"
    assert response_status(None) == "unknown"


def test_profile_requires_explicit_boolean_declarations(tmp_path):
    profile = {
        "username": "example",
        "user_alias": "",
        "participant_name": "",
        "discord_username": "",
        "email": "",
        "affiliation": "",
        "model_tag": "",
        "anon_checkbox": False,
        "paper_checkbox": True,
        "proprietary_data_checkbox": False,
        "open_code_checkbox": False,
    }
    p = tmp_path / "profile.json"
    p.write_text(json.dumps(profile))
    assert load_profile(p)["model_tag"] == ""
    profile["paper_checkbox"] = "true"
    p.write_text(json.dumps(profile))
    with pytest.raises(ValueError):
        load_profile(p)
    profile["paper_checkbox"] = True
    profile["open_code_checkbox"] = True
    p.write_text(json.dumps(profile))
    with pytest.raises(ValueError):
        load_profile(p)
