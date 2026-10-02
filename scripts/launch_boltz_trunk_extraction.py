"""Launch a recorded boltz-trunk extraction from a prepared input run."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

ACCESSIONS = {
    "cyp1a2": "P05177",
    "cyp2c9": "P11712",
    "cyp2d6": "P10635",
    "cyp3a4": "P08684",
}
REFERENCE_RUNS = {
    "cyp1a2": "20260928T111653134056Z-84ffca6a",
    "cyp2c9": "20260928T111758347957Z-0357cb7c",
    "cyp2d6": "20260928T111758358079Z-0d42000d",
    "cyp3a4": "20260928T111713518371Z-c156bd1a",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", choices=sorted(ACCESSIONS))
    parser.add_argument("prepared_run")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--recycling-steps", type=int, choices=(1, 3), default=3)
    parser.add_argument("--start-record", type=int, default=1)
    parser.add_argument("--end-record", type=int)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/shared"
    prepared = runs / args.prepared_run / "outputs/inputs"
    records = sorted(prepared.glob("*.yaml"))
    if not records:
        raise ValueError(f"No prepared YAML records found in {prepared}")
    numbered_records = [
        (int(record.name.split("_", 1)[0]), record)
        for record in records
    ]
    end_record = args.end_record or numbered_records[-1][0]
    if args.start_record < 1 or end_record < args.start_record:
        raise ValueError(
            f"Invalid record range: {args.start_record}..{end_record}"
        )
    selected_records = [
        (number, record)
        for number, record in numbered_records
        if args.start_record <= number <= end_record
    ]
    expected_count = end_record - args.start_record + 1
    if len(selected_records) != expected_count:
        raise ValueError(
            f"Expected {expected_count} records in range "
            f"{args.start_record}..{end_record}, found {len(selected_records)}"
        )
    reference = runs / REFERENCE_RUNS[args.target] / "outputs"
    if args.recycling_steps == 1:
        config_name = "extract-smoke-r1-v1.json" if args.smoke else "extract-r1-v1.json"
    else:
        config_name = "prepare-smoke-v1.json" if args.smoke else "prepare-v1.json"
    config = root / "configs/boltz-trunk" / config_name
    boltz_root = Path("/home/nagaet/ghq/github.com/N283T/boltz-trunk")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        (
            "Extract post-trunk float16 CYP/HEME/ligand tensors "
            f"for {args.target} records {args.start_record}..{end_record} "
            f"with recycle={args.recycling_steps}"
        ),
        "--config",
        str(config),
        "--input",
        f"target.a3m={reference / f'{ACCESSIONS[args.target]}.a3m'}",
        "--input",
        f"boltz_trunk_runtime.py={boltz_root / 'src/boltz/trunk_runtime.py'}",
        "--input",
        f"boltz_main.py={boltz_root / 'src/boltz/main.py'}",
        "--input",
        f"boltz_model.py={boltz_root / 'src/boltz/model/models/boltz2.py'}",
        "--input",
        f"boltz_schema.py={boltz_root / 'src/boltz/data/parse/schema.py'}",
        "--input",
        f"boltz_yaml.py={boltz_root / 'src/boltz/data/parse/yaml.py'}",
        "--input",
        f"boltz_pyproject.toml={boltz_root / 'pyproject.toml'}",
        "--input",
        f"boltz_uv.lock={boltz_root / 'uv.lock'}",
    ]
    for number, record in selected_records:
        command.extend(["--input", f"record_{number:06d}.yaml={record}"])
    command.extend(["--", "python", "scripts/run_boltz_trunk_extraction.py"])
    subprocess.run(command, cwd=root, check=True)


if __name__ == "__main__":
    main()
