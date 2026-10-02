"""Durable supervisor: finish four pooling runs, six readout screens, and report."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"]


def run_id(log):
    for line in log.read_text().splitlines():
        if line.startswith("run_id="):
            return line.split("=", 1)[1]
    return None


def await_run(identifier):
    path = ROOT / "artifacts/runs/shared" / identifier / "result.json"
    while not path.exists():
        time.sleep(5)
    result = json.loads(path.read_text())
    if result["status"] != "completed":
        raise RuntimeError(f"Pool failed: {identifier}: {result}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing-pool", action="append", default=[])
    args = parser.parse_args()
    existing = dict(item.split("=", 1) for item in args.existing_pool)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    directory = ROOT / "artifacts/boltz-trunk/experiments" / stamp
    directory.mkdir(parents=True)
    (directory / "plan.json").write_text(
        json.dumps(
            {
                "targets": TARGETS,
                "existing": existing,
                "regression_fits": 4 * 25 * 3 * 5,
                "classification_fits": 2 * 25 * 2 * 5,
                "folds": "fixed UMAP 5 folds",
                "tabpfn": "2.6, 8 ensemble members",
                "upload": False,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"suite_dir={directory}", flush=True)
    processes = {}
    streams = []
    pools, screens, failures = dict(existing), [], []
    for target in TARGETS:
        if target in existing:
            continue
        log = directory / f"pool-{target}.log"
        stream = log.open("x")
        streams.append(stream)
        processes[target] = (
            subprocess.Popen(
                ["pixi", "run", "python", "scripts/launch_cyp_boltz_pooling.py", target],
                cwd=ROOT,
                stdout=stream,
                stderr=subprocess.STDOUT,
            ),
            log,
        )
    for target in TARGETS:
        try:
            if target in processes:
                process, log = processes[target]
                code = process.wait()
                identifier = run_id(log)
                if code or identifier is None:
                    raise RuntimeError(f"Pool launcher failed: {target}; {log}")
                pools[target] = identifier
            await_run(pools[target])
            (directory / "pool-runs.json").write_text(json.dumps(pools, indent=2) + "\n")
            print(f"pool_completed={target}:{pools[target]}", flush=True)
            tracks = (
                ["regression", "classification"]
                if target in ["cyp2d6", "cyp3a4"]
                else ["regression"]
            )
            for track in tracks:
                log = directory / f"screen-{target}-{track}.log"
                with log.open("x") as stream:
                    code = subprocess.run(
                        [
                            "pixi",
                            "run",
                            "python",
                            "scripts/launch_cyp_boltz_screen.py",
                            target,
                            pools[target],
                            "--track",
                            track,
                        ],
                        cwd=ROOT,
                        stdout=stream,
                        stderr=subprocess.STDOUT,
                    ).returncode
                identifier = run_id(log)
                screens.append(
                    {"target": target, "track": track, "run": identifier, "exit_code": code}
                )
                (directory / "screen-runs.json").write_text(json.dumps(screens, indent=2) + "\n")
                if code:
                    failures.append({"target": target, "track": track, "log": str(log)})
                print(f"screen_finished={target}:{track}:{identifier}:exit={code}", flush=True)
        except Exception as exc:
            failures.append({"target": target, "error": repr(exc)})
            print(f"target_failed={target}:{exc}", flush=True)
    for stream in streams:
        stream.close()
    (directory / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    subprocess.run(
        ["pixi", "run", "python", "scripts/report_cyp_boltz_suite.py", str(directory)],
        cwd=ROOT,
        check=True,
    )
    if failures:
        raise RuntimeError(f"Suite completed with failures: {directory}")
    print(f"suite_completed={directory}", flush=True)


if __name__ == "__main__":
    main()
