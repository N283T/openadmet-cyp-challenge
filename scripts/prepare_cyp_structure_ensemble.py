"""Build an AFDB-aligned CYP/HEME/ligand ensemble and pocket definition."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import polars as pl


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_sequence(path: Path) -> str:
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Target sequence input is not FASTA")
    return "".join(lines[1:]).upper()


def heavy_atoms(residue: gemmi.Residue) -> list[gemmi.Atom]:
    return [atom for atom in residue if atom.element.name != "H"]


def minimum_distance(left: list[gemmi.Atom], right: list[gemmi.Atom]) -> float:
    if not left or not right:
        return math.inf
    return min(first.pos.dist(second.pos) for first in left for second in right)


def reference_mapping(
    sequence: str, polymer: gemmi.ResidueSpan
) -> tuple[gemmi.AlignmentResult, dict[int, int]]:
    alignment = gemmi.align_sequence_to_polymer(
        gemmi.expand_one_letter_sequence(sequence, gemmi.ResidueKind.AA),
        polymer,
        polymer.check_polymer_type(),
    )
    query_index = 0
    target_index = 0
    mapping = {}
    for count_text, operation in re.findall(r"(\d+)([MID])", alignment.cigar_str()):
        count = int(count_text)
        if operation == "M":
            for offset in range(count):
                mapping[target_index + offset] = query_index + offset + 1
            query_index += count
            target_index += count
        elif operation == "I":
            query_index += count
        else:
            target_index += count
    if target_index != len(polymer) or query_index != len(sequence):
        raise ValueError("Unexpected sequence-alignment CIGAR coverage")
    return alignment, mapping


def clone_polymer_chain(chain: gemmi.Chain) -> gemmi.Chain:
    cloned = gemmi.Chain(chain.name)
    for residue in chain.get_polymer():
        cloned.add_residue(residue)
    return cloned


def append_residue(chain: gemmi.Chain, residue: gemmi.Residue, sequence_number: int) -> None:
    """Append a copied non-polymer with an output-chain-unique author residue number."""
    chain.add_residue(residue)
    chain[-1].seqid = gemmi.SeqId(sequence_number, " ")


def block_from_chain(name: str, chain: gemmi.Chain) -> gemmi.cif.Block:
    structure = gemmi.Structure()
    structure.name = name
    model = structure.add_model(gemmi.Model("1"))
    model.add_chain(chain)
    structure.setup_entities()
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
    return structure.make_mmcif_document(groups).sole_block()


def ligand_reason(
    residue: gemmi.Residue,
    protein_atoms: list[gemmi.Atom],
    heme_iron_atoms: list[gemmi.Atom],
    config: dict,
) -> tuple[str, int, float, float]:
    atoms = heavy_atoms(residue)
    protein_distance = minimum_distance(atoms, protein_atoms)
    heme_distance = minimum_distance(atoms, heme_iron_atoms)
    if residue.entity_type == gemmi.EntityType.Water or residue.name in {"HOH", "DOD"}:
        reason = "water"
    elif residue.name in set(config["artifact_ccd_ids"]):
        reason = "configured_artifact_or_cofactor"
    elif len(atoms) < config["min_ligand_heavy_atoms"]:
        reason = "too_few_heavy_atoms"
    elif protein_distance > config["protein_ligand_contact_angstrom"]:
        reason = "no_protein_contact"
    elif heme_distance > config["max_ligand_heme_distance_angstrom"]:
        reason = "outside_heme_active_site"
    else:
        reason = "accepted_active_site_ligand"
    return reason, len(atoms), protein_distance, heme_distance


def safe_block_name(target: str, pdb_id: str, chain_name: str, ordinal: int) -> str:
    raw = f"{target}_{pdb_id}_{chain_name or 'blank'}_{ordinal}"
    return re.sub(r"[^A-Za-z0-9_.-]", "_", raw)


def main() -> None:  # noqa: C901, PLR0912, PLR0915
    target = json.loads((Path(os.environ["CYP_INPUT_DIR"]) / "target_config.json").read_text())
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    sequence = read_sequence(input_dir / "target.fasta")
    if len(sequence) != target["expected_length"]:
        raise ValueError("Target FASTA length differs from target config")
    sources = pl.read_csv(input_dir / "source_manifest.csv").sort("pdb_id")

    reference = gemmi.read_structure(str(input_dir / "afdb_model.cif"))
    if len(reference) != 1 or len(reference[0]) != 1:
        raise ValueError("Expected one AFDB model and chain")
    reference_polymer = reference[0][0].get_polymer()
    reference_sequence = gemmi.one_letter_code([residue.name for residue in reference_polymer])
    if reference_sequence != sequence:
        raise ValueError("AFDB coordinates and target FASTA differ")
    protein_type = reference_polymer.check_polymer_type()

    document = gemmi.cif.Document()
    reference_chain = clone_polymer_chain(reference[0][0])
    reference_chain.name = "A"
    document.add_copied_block(block_from_chain(f"{target['target']}_AFDB", reference_chain))
    block_rows = [
        {
            "target": target["target"],
            "pdb_id": "AFDB",
            "auth_chain": "A",
            "data_block": f"{target['target']}_AFDB",
            "polymer_residues": len(reference_chain),
            "atoms": sum(len(residue) for residue in reference_chain),
            "matched_ca": len(reference_chain),
            "sequence_identity_percent": 100.0,
            "core_rmsd_angstrom": 0.0,
            "heme_ccd_id": None,
            "heme_source_chain": None,
            "retained_ligand_ccd_ids": None,
            "retained_ligands": 0,
            "source_sha256": sha256(input_dir / "afdb_model.cif"),
        }
    ]
    ligand_rows = []
    contact_rows = []
    holo_chain_ids: set[str] = set()
    contact_chains: dict[int, set[str]] = defaultdict(set)
    contact_pdbs: dict[int, set[str]] = defaultdict(set)
    residue_names: dict[int, Counter[str]] = defaultdict(Counter)

    for source in sources.iter_rows(named=True):
        pdb_id = source["pdb_id"]
        source_path = input_dir / f"pdb_{pdb_id}.cif"
        if sha256(source_path) != source["sha256"]:
            raise ValueError(f"Source hash mismatch for {pdb_id}")
        structure = gemmi.read_structure(str(source_path))
        structure.setup_entities()
        if len(structure) != 1:
            raise ValueError(f"Expected one coordinate model in {pdb_id}")
        model = structure[0]
        heme_residues = [
            (chain.name, residue)
            for chain in model
            for residue in chain
            if residue.name in set(config["heme_ccd_ids"])
        ]
        heme_iron_atoms = [
            atom
            for _, residue in heme_residues
            for atom in residue
            if atom.element.name == "Fe"
        ]
        if not heme_iron_atoms:
            raise ValueError(f"Selected CYP structure lacks HEM iron: {pdb_id}")
        nonpolymers = [
            (chain.name, residue)
            for chain in model
            for residue in chain
            if residue.entity_type in {gemmi.EntityType.NonPolymer, gemmi.EntityType.Water}
            or residue.het_flag != "A"
        ]
        accepted_chain_count = 0
        for chain in model:
            polymer = chain.get_polymer()
            if len(polymer) < config["min_matched_ca"]:
                continue
            alignment, mapping = reference_mapping(sequence, polymer)
            identity = alignment.calculate_identity()
            if (
                alignment.match_count < config["min_matched_ca"]
                or identity < config["min_sequence_identity_percent"]
            ):
                continue
            superposition = gemmi.calculate_superposition(
                reference_polymer,
                polymer,
                protein_type,
                gemmi.SupSelect.CaP,
                trim_cycles=config["trim_cycles"],
                trim_cutoff=config["trim_cutoff"],
            )
            core_rmsd = superposition.rmsd * math.sqrt(3)
            if (
                superposition.count < config["min_matched_ca"]
                or core_rmsd > config["max_core_rmsd_angstrom"]
            ):
                continue
            accepted_chain_count += 1
            chain_id = f"{pdb_id}:{chain.name}:{accepted_chain_count}"
            protein_atoms = [atom for residue in polymer for atom in heavy_atoms(residue)]
            same_chain_hemes = [item for item in heme_residues if item[0] == chain.name]
            heme_candidates = same_chain_hemes or heme_residues
            heme_chain, selected_heme = min(
                heme_candidates,
                key=lambda item: minimum_distance(heavy_atoms(item[1]), protein_atoms),
            )
            accepted_ligands = []
            for ligand_chain, residue in nonpolymers:
                reason, atom_count, protein_distance, heme_distance = ligand_reason(
                    residue, protein_atoms, heme_iron_atoms, config
                )
                if reason != "water":
                    ligand_rows.append(
                        {
                            "target": target["target"],
                            "pdb_id": pdb_id,
                            "protein_chain": chain.name,
                            "ligand_chain": ligand_chain,
                            "ligand_seqid": str(residue.seqid),
                            "ccd_id": residue.name,
                            "heavy_atoms": atom_count,
                            "min_protein_distance_angstrom": protein_distance,
                            "min_heme_iron_distance_angstrom": heme_distance,
                            "decision": reason,
                        }
                    )
                if reason == "accepted_active_site_ligand":
                    accepted_ligands.append((ligand_chain, residue))
            if accepted_ligands:
                holo_chain_ids.add(chain_id)
            for polymer_index, residue in enumerate(polymer):
                reference_position = mapping.get(polymer_index)
                if reference_position is None:
                    continue
                residue_atoms = heavy_atoms(residue)
                for ligand_chain, ligand in accepted_ligands:
                    distance = minimum_distance(residue_atoms, heavy_atoms(ligand))
                    if distance <= config["protein_ligand_contact_angstrom"]:
                        contact_chains[reference_position].add(chain_id)
                        contact_pdbs[reference_position].add(pdb_id)
                        residue_names[reference_position][residue.name] += 1
                        contact_rows.append(
                            {
                                "target": target["target"],
                                "pdb_id": pdb_id,
                                "protein_chain": chain.name,
                                "ligand_chain": ligand_chain,
                                "ligand_seqid": str(ligand.seqid),
                                "ligand_ccd_id": ligand.name,
                                "uniprot_position": reference_position,
                                "pdb_residue_name": residue.name,
                                "pdb_residue_seqid": str(residue.seqid),
                                "min_distance_angstrom": distance,
                            }
                        )

            cloned = clone_polymer_chain(chain)
            next_seqid = max(residue.seqid.num for residue in cloned) + 1
            append_residue(cloned, selected_heme, next_seqid)
            next_seqid += 1
            for _, ligand in accepted_ligands:
                append_residue(cloned, ligand, next_seqid)
                next_seqid += 1
            fragment = gemmi.Structure()
            fragment.name = "fragment"
            fragment_model = fragment.add_model(gemmi.Model("1"))
            moved = fragment_model.add_chain(cloned)
            fragment_model.transform_pos_and_adp(superposition.transform)
            moved.name = chain.name or "_"
            block_name = safe_block_name(
                target["target"], pdb_id.upper(), chain.name, accepted_chain_count
            )
            document.add_copied_block(block_from_chain(block_name, moved))
            block_rows.append(
                {
                    "target": target["target"],
                    "pdb_id": pdb_id,
                    "auth_chain": chain.name,
                    "data_block": block_name,
                    "polymer_residues": len(polymer),
                    "atoms": sum(len(residue) for residue in cloned),
                    "matched_ca": superposition.count,
                    "sequence_identity_percent": identity,
                    "core_rmsd_angstrom": core_rmsd,
                    "heme_ccd_id": selected_heme.name,
                    "heme_source_chain": heme_chain,
                    "retained_ligand_ccd_ids": ";".join(
                        ligand.name for _, ligand in accepted_ligands
                    ),
                    "retained_ligands": len(accepted_ligands),
                    "source_sha256": source["sha256"],
                }
            )
        if accepted_chain_count == 0:
            raise ValueError(f"No target-like protein chain found in {pdb_id}")

    ensemble_name = f"{target['target'].lower()}_aligned_complex_ensemble.cif"
    ensemble_path = output_dir / ensemble_name
    document.write_file(str(ensemble_path))
    checked = gemmi.cif.read_file(str(ensemble_path))
    if len(checked) != len(block_rows):
        raise ValueError("Multi-block ensemble count changed during write")
    for block, row in zip(checked, block_rows, strict=True):
        parsed = gemmi.make_structure_from_block(block)
        if len(parsed) != 1 or len(parsed[0]) != 1:
            raise ValueError(f"Block is not a single-chain structure: {block.name}")
        parsed_chain = parsed[0][0]
        retained_nonpolymers = [
            residue
            for residue in parsed_chain
            if residue.entity_type != gemmi.EntityType.Polymer
        ]
        retained_hemes = [
            residue
            for residue in retained_nonpolymers
            if residue.name in set(config["heme_ccd_ids"])
        ]
        retained_ligands = [
            residue
            for residue in retained_nonpolymers
            if residue.name not in set(config["heme_ccd_ids"])
        ]
        expected_hemes = 0 if row["pdb_id"] == "AFDB" else 1
        if len(retained_hemes) != expected_hemes:
            raise ValueError(f"HEME count changed in block: {block.name}")
        if len(retained_ligands) != row["retained_ligands"]:
            raise ValueError(f"Ligand count changed in block: {block.name}")
        if parsed[0].count_atom_sites() != row["atoms"]:
            raise ValueError(f"Atom count changed in block: {block.name}")

    holo_count = len(holo_chain_ids)
    if holo_count == 0:
        raise ValueError("No active-site ligands survived artifact filtering")
    pocket_rows = []
    core_positions = []
    extended_positions = []
    for position in sorted(contact_chains):
        chain_count = len(contact_chains[position])
        fraction = chain_count / holo_count
        tier = "all"
        is_core = fraction >= config["core_contact_fraction"]
        is_extended = (
            is_core
            or fraction >= config["extended_contact_fraction"]
            and chain_count >= config["extended_min_holo_chains"]
        )
        if is_core:
            tier = "core"
            core_positions.append(position)
        if is_extended:
            extended_positions.append(position)
            if tier == "all":
                tier = "extended"
        pocket_rows.append(
            {
                "target": target["target"],
                "accession": target["accession"],
                "uniprot_position": position,
                "reference_amino_acid": sequence[position - 1],
                "contact_holo_chains": chain_count,
                "contact_pdbs": len(contact_pdbs[position]),
                "holo_chain_fraction": fraction,
                "observed_pdb_residue_names": ";".join(sorted(residue_names[position])),
                "tier": tier,
            }
        )

    pl.DataFrame(block_rows).write_csv(output_dir / "block_map.csv")
    pl.DataFrame(ligand_rows).write_csv(output_dir / "ligand_inventory.csv")
    pl.DataFrame(contact_rows).write_csv(output_dir / "pocket_contacts.csv")
    pl.DataFrame(pocket_rows).write_csv(output_dir / "pocket_residues.csv")
    pocket_definition = {
        "target": target["target"],
        "accession": target["accession"],
        "numbering": "1-based canonical UniProt/AFDB sequence positions",
        "contact_cutoff_angstrom": config["protein_ligand_contact_angstrom"],
        "ligand_rule": (
            "non-water, non-HEM, non-configured-artifact residue with at least "
            f"{config['min_ligand_heavy_atoms']} heavy atoms, within "
            f"{config['max_ligand_heme_distance_angstrom']} A of HEM iron"
        ),
        "holo_chains": holo_count,
        "core_min_fraction": config["core_contact_fraction"],
        "core_residues": core_positions,
        "extended_min_fraction": config["extended_contact_fraction"],
        "extended_min_holo_chains": config["extended_min_holo_chains"],
        "extended_residues": extended_positions,
        "all_contact_residues": sorted(contact_chains),
        "artifact_ccd_ids": config["artifact_ccd_ids"],
    }
    (output_dir / "pocket_residues.json").write_text(
        json.dumps(pocket_definition, indent=2) + "\n"
    )
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "target": target["target"],
        "accession": target["accession"],
        "source_structures": sources.height,
        "experimental_chains": len(block_rows) - 1,
        "data_blocks": len(block_rows),
        "protein_chains_per_block": 1,
        "afdb_block_has_heme": False,
        "retained_heme_residues": len(block_rows) - 1,
        "retained_active_site_ligands": sum(row["retained_ligands"] for row in block_rows),
        "retained_water_or_configured_artifact_residues": 0,
        "holo_chains_for_pocket": holo_count,
        "accepted_ligand_instances": sum(
            row["decision"] == "accepted_active_site_ligand" for row in ligand_rows
        ),
        "core_pocket_residues": len(core_positions),
        "extended_pocket_residues": len(extended_positions),
        "all_contact_residues": len(contact_chains),
        "ensemble_file": ensemble_name,
        "ensemble_bytes": ensemble_path.stat().st_size,
        "ensemble_sha256": sha256(ensemble_path),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
