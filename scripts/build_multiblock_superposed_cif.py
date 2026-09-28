"""Write one mmCIF with AFDB and each fitted chain in its own data_ block."""

from __future__ import annotations

import hashlib
import json
import math
import os
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import polars as pl


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residual_transform(transform: gemmi.Transform) -> tuple[float, float]:
    translation = math.sqrt(sum(value * value for value in transform.vec.tolist()))
    rotation = max(
        abs(value - (1.0 if row == column else 0.0))
        for row, values in enumerate(transform.mat.tolist())
        for column, value in enumerate(values)
    )
    return translation, rotation


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    source_path = input_dir / "all_chains.cif"
    source_summary = json.loads((input_dir / "all_chains_summary.json").read_text())
    if sha256(source_path) != source_summary["sha256"]:
        raise ValueError("All-chain source CIF hash mismatch")
    source = gemmi.read_structure(str(source_path))
    if len(source) != 1 or len(source[0]) != config["expected_total_blocks"]:
        raise ValueError("Unexpected source model/chain count")
    source.setup_entities()
    reference_polymer = source[0]["AFDB_A"].get_polymer()
    if len(reference_polymer) != 503:
        raise ValueError("AFDB reference is not full-length")
    protein_type = reference_polymer.check_polymer_type()
    alignments = pl.read_csv(input_dir / "chain_alignment.csv")
    if alignments.height != config["expected_source_chains"]:
        raise ValueError("Unexpected source chain count")

    entries = [
        {"pdb_id": "AFDB", "original_chain": "A", "overlay_chain": "AFDB_A"},
        *alignments.sort("pdb_id", "original_chain").iter_rows(named=True),
    ]
    groups = gemmi.MmcifOutputGroups(
        True,
        cell=False,
        symmetry=False,
        scale=False,
        origx=False,
        assembly=False,
        ncs=False,
        diffrn=False,
        refine=False,
        reflns=False,
    )
    document = gemmi.cif.Document()
    rows = []
    for entry in entries:
        block_name = entry["overlay_chain"]
        original_chain = entry["original_chain"]
        fragment = gemmi.Structure()
        fragment.name = block_name
        model = fragment.add_model(gemmi.Model("1"))
        chain = model.add_chain(source[0][block_name])
        chain.name = original_chain
        for residue in chain:
            prefix = f"{block_name}_"
            if residue.subchain.startswith(prefix):
                residue.subchain = residue.subchain[len(prefix) :]
            if residue.entity_id.startswith(prefix):
                residue.entity_id = residue.entity_id[len(prefix) :]
        fragment.setup_entities()
        block = fragment.make_mmcif_document(groups).sole_block()
        if block.name != block_name:
            raise ValueError(f"Unexpected data block name: {block.name}")
        document.add_copied_block(block)
        rows.append(
            {
                "pdb_id": entry["pdb_id"],
                "original_chain": original_chain,
                "data_block": block_name,
                "atoms": model.count_atom_sites(),
                "heme_residues": sum(residue.name == "HEM" for residue in chain),
            }
        )

    output_path = output_dir / "cyp3a4_all_chains_143_data_blocks.cif"
    document.write_file(str(output_path))
    checked = gemmi.cif.read_file(str(output_path))
    if len(checked) != config["expected_total_blocks"]:
        raise ValueError("Output does not contain 143 data blocks")
    for block, row in zip(checked, rows, strict=True):
        if block.name != row["data_block"]:
            raise ValueError("Output data block order/name mismatch")
        structure = gemmi.make_structure_from_block(block)
        if len(structure) != 1 or len(structure[0]) != 1:
            raise ValueError(f"Data block does not contain one chain: {block.name}")
        if structure[0][0].name != row["original_chain"]:
            raise ValueError(f"Original chain ID not restored in {block.name}")
        if structure[0].count_atom_sites() != row["atoms"]:
            raise ValueError(f"Atom count changed in {block.name}")
        hemes = sum(residue.name == "HEM" for residue in structure[0][0])
        if hemes != row["heme_residues"]:
            raise ValueError(f"HEME count changed in {block.name}")
        structure.setup_entities()
        if block.name != "AFDB_A":
            fit = gemmi.calculate_superposition(
                reference_polymer,
                structure[0][0].get_polymer(),
                protein_type,
                gemmi.SupSelect.CaP,
                trim_cycles=5,
                trim_cutoff=2,
            )
            translation, rotation = residual_transform(fit.transform)
            if (
                fit.count < config["min_matched_ca"]
                or fit.rmsd * math.sqrt(3) > config["max_core_rmsd_angstrom"]
                or translation > config["max_residual_translation_angstrom"]
                or rotation > config["max_residual_rotation_element"]
            ):
                raise ValueError(f"Chain not independently fitted to AFDB: {block.name}")
            row["matched_ca"] = fit.count
            row["core_rmsd_angstrom"] = fit.rmsd * math.sqrt(3)
            row["residual_translation_angstrom"] = translation
            row["residual_rotation_element"] = rotation
        else:
            if len(structure[0][0].get_polymer()) != 503:
                raise ValueError("AFDB block is not full-length")
            row["matched_ca"] = 503
            row["core_rmsd_angstrom"] = 0.0
            row["residual_translation_angstrom"] = 0.0
            row["residual_rotation_element"] = 0.0

    pl.DataFrame(rows).write_csv(output_dir / "block_map.csv")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "all_chain_run_id": config["all_chain_run_id"],
        "data_blocks": len(checked),
        "source_chains_fitted": len(checked) - 1,
        "chains_per_block": 1,
        "atoms_across_blocks": sum(row["atoms"] for row in rows),
        "heme_residues_across_blocks": sum(row["heme_residues"] for row in rows),
        "max_residual_translation_angstrom": max(
            row["residual_translation_angstrom"] for row in rows
        ),
        "max_residual_rotation_element": max(row["residual_rotation_element"] for row in rows),
        "file": output_path.name,
        "bytes": output_path.stat().st_size,
        "sha256": sha256(output_path),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
