"""Record a report and promote manageable diagnostic outputs with hashes."""

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_id")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = root / "artifacts/runs/classification" / args.run_id
    if json.loads((source / "result.json").read_text())["status"] != "completed":
        raise ValueError("Only completed diagnostics can be retained")
    config = root / "configs/boltz-trunk" / f"tdi-diagnostics-report-{args.run_id}.json"
    with config.open("x") as stream:
        json.dump({"source_run": args.run_id}, stream, indent=2)
        stream.write("\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Render labelled-only TDI trunk integration diagnostics and limitations",
        "--config",
        str(config),
    ]
    for path in sorted((source / "outputs").glob("*")):
        if path.is_file():
            command.extend(["--input", f"{path.name}={path}"])
    result = subprocess.run(
        command + ["--", "python", "scripts/report_cyp_boltz_tdi_diagnostics.py"],
        cwd=root,
        text=True,
        capture_output=True,
    )
    print(result.stdout)
    if result.returncode:
        raise RuntimeError(result.stderr)
    report_id = next(
        line.split("=", 1)[1] for line in result.stdout.splitlines() if line.startswith("run_id=")
    )
    report = root / "artifacts/runs/classification" / report_id
    tracked = root / "artifacts/tracked" / f"cyp-boltz-tdi-diagnostics-{args.run_id}"
    tracked.mkdir()
    for path in sorted((source / "outputs").glob("*")):
        if path.is_file():
            shutil.copy2(path, tracked / path.name)
    shutil.copy2(report / "outputs/report.md", tracked / "report.md")
    inventories = {
        path.name: {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
        for path in sorted(tracked.iterdir())
        if path.is_file()
    }
    (tracked / "provenance.json").write_text(
        json.dumps(
            {
                "diagnostic_run": args.run_id,
                "report_run": report_id,
                "source_inputs": json.loads((source / "manifest.json").read_text())["inputs"],
                "files": inventories,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"tracked={tracked}")
    print(f"report_run={report_id}")


if __name__ == "__main__":
    main()
