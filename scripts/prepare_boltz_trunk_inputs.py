"""Expand one fixed CYP/HEME context and a compound table into Boltz inputs."""

from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path

import polars as pl
import yaml


def sha256(path: Path) -> str:
    """Return a file checksum for immutable input provenance."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> str:
    """Read one FASTA sequence."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("target.fasta must contain one FASTA record")
    sequence = "".join(lines[1:]).replace(" ", "").upper()
    if not sequence or not re.fullmatch(r"[A-Z]+", sequence):
        raise ValueError("target.fasta contains an invalid protein sequence")
    return sequence


def read_a3m_query(path: Path) -> str:
    """Read and normalize the first A3M sequence."""
    sequence_parts: list[str] = []
    seen_header = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header and sequence_parts:
                    break
                seen_header = True
                continue
            if seen_header:
                sequence_parts.append(line)
    query = "".join(sequence_parts)
    return "".join(char for char in query if not char.islower() and char not in "-.").upper()


def safe_record_id(compound_id: str, index: int) -> str:
    """Make a deterministic Boltz record/file identifier."""
    suffix = re.sub(r"[^A-Za-z0-9_.-]+", "_", compound_id).strip("._-") or "compound"
    return f"{index:06d}_{suffix[:80]}"


def validate_config(config: dict[str, object]) -> None:
    """Reject settings that would silently switch to another modeling route."""
    chain_ids = [
        config["protein_chain_id"],
        config["heme_chain_id"],
        config["ligand_chain_id"],
    ]
    if len(set(chain_ids)) != 3:
        raise ValueError("protein, HEM, and ligand chain IDs must be distinct")
    if config["heme_ccd_id"] != "HEM":
        raise ValueError("The v1 fixed cofactor must use CCD HEM")
    if config["templates"] != "none":
        raise ValueError("The primary v1 preparation is template-free")
    if config["checkpoint"] != "boltz2_conf.ckpt":
        raise ValueError("Trunk extraction must use the structure checkpoint")
    if config["device"] != "cuda" or config["cuda_kernels"] is not True:
        raise ValueError("The v1 runtime requires CUDA with kernels enabled")
    if config["compound_scope"] != "target_labeled_train_plus_test":
        raise ValueError("The v1 input scope must be target-labeled train plus test")
    if config["embedding_stage"] != "post":
        raise ValueError("The v1 extraction must retain post-trunk tensors only")
    if config["output_dtype"] != "float16":
        raise ValueError("The v1 extraction must store tensors as float16")


def select_target_compounds(
    input_dir: Path,
    target: str,
) -> pl.DataFrame:
    """Select target-labeled train compounds and every valid test compound."""
    compounds = pl.read_parquet(input_dir / "compounds.parquet").select(
        "compound_id", "standardized_smiles"
    )
    aliases = pl.read_parquet(input_dir / "aliases.parquet").select(
        "record_id", "role", "compound_id"
    )
    measurement_sources = pl.read_parquet(
        input_dir / "measurement_sources.parquet"
    ).select("measurement_id", "record_id")
    measurements = pl.read_parquet(input_dir / "measurements.parquet").select(
        "measurement_id", "compound_id", "endpoint"
    )
    target_train = (
        measurements.join(measurement_sources, on="measurement_id", how="inner")
        .join(
            aliases.select("record_id", "role").unique(),
            on="record_id",
            how="inner",
        )
        .filter(
            (pl.col("role") == "train")
            & pl.col("endpoint").str.starts_with(f"{target}_")
        )
        .select("compound_id")
        .unique()
        .with_columns(pl.lit("train").alias("role"))
    )
    test = (
        aliases.filter(pl.col("role") == "test")
        .select("compound_id")
        .unique()
        .with_columns(pl.lit("test").alias("role"))
    )
    selected = pl.concat([target_train, test]).filter(pl.col("compound_id").is_not_null())
    if selected["compound_id"].n_unique() != selected.height:
        raise ValueError("Target train and test compound sets overlap")
    result = selected.join(compounds, on="compound_id", how="inner", validate="1:1").sort(
        "role", "compound_id"
    )
    if result.height != selected.height:
        raise ValueError("Selected train/test compounds are missing standardized SMILES")
    return result


def main() -> None:
    """Write standard Boltz JSON records and an auditable manifest."""
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    config_path = Path(os.environ["CYP_CONFIG_PATH"])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    validate_config(config)

    fasta_path = input_dir / "target.fasta"
    msa_path = input_dir / "target.a3m"
    pocket_path = input_dir / "pocket_residues.json"
    sequence = read_fasta(fasta_path)
    if read_a3m_query(msa_path) != sequence:
        raise ValueError("MSA query does not match target FASTA")

    pocket = json.loads(pocket_path.read_text(encoding="utf-8"))
    residue_set = str(config["pocket_residue_set"])
    residues = pocket.get(residue_set)
    if not isinstance(residues, list) or not residues or not all(
        isinstance(value, int) and not isinstance(value, bool) for value in residues
    ):
        raise ValueError(f"Pocket file has no valid {residue_set}")
    if min(residues) < 1 or max(residues) > len(sequence):
        raise ValueError("Pocket residues fall outside the canonical target sequence")

    compounds = select_target_compounds(input_dir, str(pocket["target"]))
    if compounds["compound_id"].n_unique() != compounds.height:
        raise ValueError("Compound IDs must be unique")
    if compounds.null_count().sum_horizontal().item() != 0:
        raise ValueError("Compound table contains null IDs or SMILES")
    limit = config.get("limit")
    if limit is not None:
        limit = int(limit)
        train_limit = (limit + 1) // 2
        test_limit = limit // 2
        compounds = pl.concat(
            [
                compounds.filter(pl.col("role") == "train").head(train_limit),
                compounds.filter(pl.col("role") == "test").head(test_limit),
            ]
        )

    records_dir = output_dir / "inputs"
    records_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    protein_id = str(config["protein_chain_id"])
    heme_id = str(config["heme_chain_id"])
    ligand_id = str(config["ligand_chain_id"])
    msa_resolved = str(msa_path.resolve())
    for index, row in enumerate(compounds.iter_rows(named=True), start=1):
        compound_id = str(row["compound_id"])
        record_id = safe_record_id(compound_id, index)
        payload = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": protein_id,
                        "sequence": sequence,
                        "msa": msa_resolved,
                    }
                },
                {"ligand": {"id": heme_id, "ccd": str(config["heme_ccd_id"])}},
                {"ligand": {"id": ligand_id, "smiles": row["standardized_smiles"]}},
            ],
            "constraints": [
                {
                    "pocket": {
                        "binder": ligand_id,
                        "contacts": [[protein_id, residue] for residue in residues],
                        "max_distance": float(config["pocket_max_distance_angstrom"]),
                        "force": bool(config["pocket_force"]),
                    }
                }
            ],
        }
        path = records_dir / f"{record_id}.yaml"
        path.write_text(
            yaml.safe_dump(payload, sort_keys=False, default_flow_style=False),
            encoding="utf-8",
        )
        rows.append(
            {
                "record_id": record_id,
                "compound_id": compound_id,
                "role": row["role"],
                "standardized_smiles": row["standardized_smiles"],
                "input_path": str(path),
                "input_sha256": sha256(path),
            }
        )

    pl.DataFrame(rows).write_parquet(output_dir / "manifest.parquet")
    runtime_args = [
        "extract-trunk",
        str(records_dir),
        "--output-dir",
        "<CYP_OUTPUT_DIR>/embeddings",
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
        runtime_args.append("--flash-attn")
    if not config["bf16"]:
        runtime_args.append("--no-bf16")
    (output_dir / "runtime_args.json").write_text(
        json.dumps(runtime_args, indent=2) + "\n", encoding="utf-8"
    )
    summary = {
        "target": pocket["target"],
        "accession": pocket["accession"],
        "records": compounds.height,
        "train_records": compounds.filter(pl.col("role") == "train").height,
        "test_records": compounds.filter(pl.col("role") == "test").height,
        "compound_scope": config["compound_scope"],
        "sequence_length": len(sequence),
        "msa_sha256": sha256(msa_path),
        "pocket_residue_set": residue_set,
        "pocket_residues": residues,
        "pocket_max_distance_angstrom": config["pocket_max_distance_angstrom"],
        "pocket_force": config["pocket_force"],
        "fixed_entities": [{"id": heme_id, "ccd": config["heme_ccd_id"]}],
        "variable_entity": {"id": ligand_id, "source": "standardized_smiles"},
        "properties": [],
        "templates": [],
        "runtime": {key: config[key] for key in (
            "checkpoint",
            "recycling_steps",
            "device",
            "cuda_kernels",
            "flash_attention",
            "bf16",
            "max_msa_sequences",
            "subsample_msa",
            "num_subsampled_msa",
            "seed",
            "ligand_conformer_seed",
            "ligand_reference_conformer",
            "embedding_stage",
            "output_dtype",
        )},
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
