"""Send prepared submissions once, with private profiles and durable uncertain outcomes."""

import json
import shutil
from pathlib import Path

from .api import SPACE, leaderboard
from .database import connect
from .storage import allocate, digest, now, snapshot_source, write_json
from .tracks import ENDPOINTS, validate

TEXT_FIELDS = {
    "username",
    "user_alias",
    "participant_name",
    "discord_username",
    "email",
    "affiliation",
    "model_tag",
}
BOOL_FIELDS = {"anon_checkbox", "paper_checkbox", "proprietary_data_checkbox", "open_code_checkbox"}
TRACK_LABELS = {
    "regression": "Regression Prediction",
    "classification": "Classification Prediction",
}


def load_profile(path):
    profile = json.loads(Path(path).read_text())
    if set(profile) != TEXT_FIELDS | BOOL_FIELDS:
        raise ValueError("Profile keys must match configs/submission-profile.example.json")
    if any(not isinstance(profile[k], str) for k in TEXT_FIELDS):
        raise ValueError("Profile text values must be strings")
    if any(type(profile[k]) is not bool for k in BOOL_FIELDS):
        raise ValueError("Profile checkbox values must be JSON booleans")
    if not profile["username"].strip():
        raise ValueError("Profile requires a Hugging Face username")
    if profile["open_code_checkbox"] and not profile["model_tag"].strip():
        raise ValueError("Public-code declaration requires a report URL")
    return profile


def response_status(response):
    value = response.get("value", "") if isinstance(response, dict) else response
    if not isinstance(value, str):
        return "unknown"
    if value.startswith("Submission received from "):
        return "sent"
    if value.startswith("Error:"):
        return "rejected"
    return "unknown"


def snapshot_leaderboards(track):
    """Save every current leaderboard before a live upload; fail closed."""

    snapshots = []
    for endpoint in ["MA", *ENDPOINTS[track]]:
        path = leaderboard(track, endpoint)
        snapshots.append({"endpoint": endpoint, "path": str(path)})
    return snapshots


def check_previous_attempts(previous, submission_id, resubmit_sent, retry_rejected):
    """Require an explicit, separate submission for a known sent payload."""

    if any(status == "unknown" for _, status in previous):
        raise ValueError("Payload delivery uncertain; inspect receipts before any resubmission")
    if any(prior_id == submission_id and status == "sent" for prior_id, status in previous):
        raise ValueError("This prepared submission has already been attempted")
    if any(status == "sent" for _, status in previous) and not resubmit_sent:
        raise ValueError("Payload already sent; use a new prepared ID and --resubmit-sent")
    if any(status == "rejected" for _, status in previous) and not retry_rejected:
        raise ValueError("Previously rejected payload; inspect receipt before --retry-rejected")


def submit(args):
    profile_path = args.profile.resolve()
    profile = load_profile(profile_path)
    with connect() as conn:
        row = conn.execute(
            "SELECT s.track,s.manifest,c.manifest FROM submissions s "
            "JOIN candidates c ON c.id=s.candidate_id WHERE s.id=%s",
            (args.id,),
        ).fetchone()
        if not row:
            raise ValueError("Unknown prepared submission")
        track, manifest, candidate = row
        file = Path(manifest["file"])
        if digest(file) != manifest["sha256"] or manifest["sha256"] != candidate["sha256"]:
            raise ValueError("Submission file changed")
        reference = Path(candidate["file"]).parent / ("reference" + Path(candidate["file"]).suffix)
        # Reference may have a different extension from the prediction file.
        references = list(Path(candidate["file"]).parent.glob("reference.*"))
        if len(references) != 1:
            raise ValueError("Missing or ambiguous reference snapshot")
        reference = references[0]
        if digest(reference) != candidate["reference_sha256"]:
            raise ValueError("Reference changed")
        report = validate(file, track, reference)
        if args.dry_run:
            print(
                json.dumps(
                    dict(
                        submission_id=args.id,
                        track=track,
                        sha256=manifest["sha256"],
                        validation=report,
                        profile_fields=sorted(profile),
                        report_url_present=bool(profile["model_tag"]),
                        sent=False,
                    )
                )
            )
            return
        # Session lock survives the commit that makes intent durable before network I/O.
        conn.execute("SELECT pg_advisory_lock(hashtext(%s))", ("cyp-submit:" + track,))
        previous = conn.execute(
            "SELECT submission_id,status FROM submission_attempts "
            "WHERE track=%s AND file_sha256=%s",
            (track, manifest["sha256"]),
        ).fetchall()
        check_previous_attempts(previous, args.id, args.resubmit_sent, args.retry_rejected)
        # These network reads happen before any upload intent or submission API call.
        # A failure propagates and prevents the live upload.
        pre_submission_leaderboards = snapshot_leaderboards(track)
        identifier, path = allocate("submission-attempts", track)
        path.chmod(0o700)
        upload = path / file.name
        shutil.copy2(file, upload)
        if digest(upload) != manifest["sha256"]:
            raise ValueError("File changed during copy")
        private = path / "profile.json"
        write_json(private, profile)
        private.chmod(0o600)
        intent = dict(
            id=identifier,
            submission_id=args.id,
            track=track,
            created_at=now(),
            file_sha256=digest(upload),
            profile_sha256=digest(private),
            endpoint="/submit_predictions",
            pre_submission_leaderboards=pre_submission_leaderboards,
            source=snapshot_source(path / "source"),
            status="unknown",
        )
        write_json(path / "intent.json", intent)
        conn.execute(
            "INSERT INTO submission_attempts "
            "(id,submission_id,track,file_sha256,status,manifest_path) "
            "VALUES (%s,%s,%s,%s,%s,%s)",
            (identifier, args.id, track, digest(upload), "unknown", str(path / "intent.json")),
        )
        conn.execute("UPDATE submissions SET status='unknown' WHERE id=%s", (args.id,))
        conn.commit()
        outcome = dict(intent)
        try:
            from gradio_client import Client, handle_file

            client = Client(SPACE, download_files=str(path / "downloads"), verbose=False)
            info = client.view_api(return_format="dict", print_info=False)
            write_json(path / "api-info.json", info)
            endpoint = info["named_endpoints"]["/submit_predictions"]
            names = {p["parameter_name"] for p in endpoint["parameters"]}
            if names != set(profile) | {"track_select", "file_input"}:
                raise ValueError("Live submission API fields changed; refusing upload")
            # Exactly one call. A timeout may follow successful delivery; never retry here.
            response = client.predict(
                **profile,
                track_select=TRACK_LABELS[track],
                file_input=handle_file(str(upload)),
                api_name="/submit_predictions",
            )
            write_json(path / "response.json", response)
            outcome.update(status=response_status(response), response_file="response.json")
        except Exception as exc:
            outcome.update(status="unknown", error_type=type(exc).__name__)
        write_json(path / "result.json", outcome)
        conn.execute(
            "UPDATE submission_attempts SET status=%s WHERE id=%s", (outcome["status"], identifier)
        )
        conn.execute("UPDATE submissions SET status=%s WHERE id=%s", (outcome["status"], args.id))
        conn.execute(
            "INSERT INTO submission_events "
            "(submission_id,status,evidence_path,evidence_sha256,notes) "
            "VALUES (%s,%s,%s,%s,%s)",
            (
                args.id,
                outcome["status"],
                str(path / "result.json"),
                digest(path / "result.json"),
                "Automated attempt " + identifier + "; sent means receipt, not completed scoring",
            ),
        )
        conn.commit()
        print(
            json.dumps(
                dict(
                    attempt=identifier, status=outcome["status"], receipt=str(path / "result.json")
                )
            )
        )
        if outcome["status"] != "sent":
            raise RuntimeError(
                "Submission not confirmed; inspect saved receipt before further action"
            )
