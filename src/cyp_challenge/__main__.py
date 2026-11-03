"""Single command entry point for infrastructure and recorded experiments."""

import argparse
from pathlib import Path

from . import api, candidates, database, runs, submission
from .tracks import ENDPOINTS


def main():
    parser = argparse.ArgumentParser(prog="cyp")
    sub = parser.add_subparsers(dest="action", required=True)
    db = sub.add_parser("db")
    db.add_argument("operation", choices=["init", "start", "stop", "migrate", "backup"])
    sub.add_parser("api-inspect")
    lb = sub.add_parser("leaderboard")
    lb.add_argument("--track", choices=ENDPOINTS, required=True)
    lb.add_argument("--endpoint", default="MA")
    lb.add_argument("--board", choices=["live", "interim"], default="live")
    sub.add_parser("register-raw")
    sub.add_parser("reconcile")
    prune = sub.add_parser("prune-runs")
    prune.add_argument("--status", choices=["failed"], default="failed")
    prune.add_argument("--execute", action="store_true")
    note = sub.add_parser("note")
    note.add_argument("id")
    note.add_argument("--decision", choices=["observation", "keep", "discard"], required=True)
    note.add_argument("--notes", required=True)
    note.add_argument("--issue")
    replay = sub.add_parser("replay")
    replay.add_argument("id")
    run = sub.add_parser("run")
    run.add_argument("--track", choices=["shared", *ENDPOINTS], required=True)
    run.add_argument("--purpose", required=True)
    run.add_argument("--issue")
    run.add_argument("--config", type=Path)
    run.add_argument(
        "--input-mode",
        choices=["copy", "hardlink"],
        default="copy",
        help=(
            "Hardlink retained immutable large inputs on the same filesystem; "
            "copies remain the default."
        ),
    )
    run.add_argument(
        "--input",
        action="append",
        default=[],
        metavar="[ALIAS=]PATH",
        help=(
            "Snapshot an input. Prefer a stable ALIAS=PATH name; bare paths retain "
            "the legacy numeric-prefix behavior."
        ),
    )
    run.add_argument("command", nargs=argparse.REMAINDER)
    candidate = sub.add_parser("candidate")
    candidate.add_argument("--track", choices=ENDPOINTS, required=True)
    candidate.add_argument("--file", type=Path, required=True)
    candidate.add_argument("--reference", type=Path, required=True)
    candidate.add_argument("--run", action="append", required=True)
    candidate.add_argument("--notes", required=True)
    prepare = sub.add_parser("prepare-submission")
    prepare.add_argument("id")
    event = sub.add_parser("submission-event")
    event.add_argument("id")
    event.add_argument(
        "--status", choices=["sent", "unknown", "accepted", "rejected"], required=True
    )
    event.add_argument("--evidence", type=Path, required=True)
    event.add_argument("--notes", required=True)
    send = sub.add_parser("submit")
    send.add_argument("id")
    send.add_argument("--profile", type=Path, default=Path(".private/submission-profile.json"))
    send.add_argument("--dry-run", action="store_true")
    send.add_argument("--retry-rejected", action="store_true")
    send.add_argument(
        "--confirm",
        help='Non-interactive confirmation "SUBMIT <id>"; use only after explicit user approval',
    )
    send.add_argument(
        "--resubmit-sent",
        action="store_true",
        help="Explicitly resubmit a payload sent under a different submission ID",
    )
    args = parser.parse_args()
    if args.action == "submit":
        submission.submit(args)
    elif args.action == "db":
        database.lifecycle(args.operation)
    elif args.action == "api-inspect":
        api.inspect_api()
    elif args.action == "leaderboard":
        api.leaderboard(args.track, args.endpoint, args.board)
    elif args.action == "register-raw":
        api.register_raw()
    elif args.action == "run":
        return runs.run(args)
    elif args.action == "note":
        runs.note(args)
    elif args.action == "replay":
        return runs.replay(args.id)
    elif args.action == "reconcile":
        runs.reconcile()
    elif args.action == "prune-runs":
        runs.prune(args)
    elif args.action == "candidate":
        candidates.candidate(args)
    elif args.action == "prepare-submission":
        candidates.prepare(args.id)
    elif args.action == "submission-event":
        candidates.event(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
