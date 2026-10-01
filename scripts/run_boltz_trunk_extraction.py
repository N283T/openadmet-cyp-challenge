"""Run the external boltz-trunk encoder inside a recorded CYP run."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import yaml


def main() -> None:
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    scratch_dir = Path(os.environ["TMPDIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text(encoding="utf-8"))
    records_dir = scratch_dir / "boltz-inputs"
    records_dir.mkdir(parents=True, exist_ok=False)

    source_records = sorted(input_dir.glob("record_*.yaml"))
    if not source_records:
        raise ValueError("No record_*.yaml inputs were supplied")
    msa_path = (input_dir / "target.a3m").resolve()
    for source in source_records:
        payload = yaml.safe_load(source.read_text(encoding="utf-8"))
        proteins = [
            entity["protein"]
            for entity in payload["sequences"]
            if "protein" in entity
        ]
        if len(proteins) != 1:
            raise ValueError(f"Expected one protein in {source}")
        proteins[0]["msa"] = str(msa_path)
        (records_dir / source.name).write_text(
            yaml.safe_dump(payload, sort_keys=False, default_flow_style=False),
            encoding="utf-8",
        )

    boltz_root = Path("/home/nagaet/ghq/github.com/N283T/boltz-trunk")
    executable = boltz_root / ".venv/bin/boltz-stfree"
    if not executable.is_file():
        raise FileNotFoundError(f"Missing boltz-trunk executable: {executable}")
    command = [
        str(executable),
        "extract-trunk",
        str(records_dir),
        "--output-dir",
        str(output_dir),
        "--recycling-steps",
        str(config["recycling_steps"]),
        "--device",
        str(config["device"]),
        "--max-msa-seqs",
        str(config["max_msa_sequences"]),
        "--subsample-msa",
        "--num-subsampled-msa",
        str(config["num_subsampled_msa"]),
        "--seed",
        str(config["seed"]),
        "--ligand-conformer-seed",
        str(config["ligand_conformer_seed"]),
        "--embedding-stage",
        str(config["embedding_stage"]),
        "--output-dtype",
        str(config["output_dtype"]),
    ]
    if config["flash_attention"]:
        command.append("--flash-attn")
    if not config["cuda_kernels"]:
        command.append("--no-kernels")
    if not config["bf16"]:
        command.append("--no-bf16")
    subprocess.run(command, cwd=boltz_root, check=True)


if __name__ == "__main__":
    main()
