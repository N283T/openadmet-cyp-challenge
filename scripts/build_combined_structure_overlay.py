"""Combine all aligned public CYP3A4 structures into one visible mmCIF model."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import polars as pl


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_structure(
    model: gemmi.Model,
    structure: gemmi.Structure,
    source_id: str,
    source_file: str,
    used_names: set[str],
) -> list[dict]:
    if len(structure) != 1:
        raise ValueError(f"Expected one coordinate model in {source_file}")
    rows = []
    for source_chain in structure[0]:
        original_name = source_chain.name or "blank"
        base_name = f"{source_id}_{original_name}"
        new_name = base_name
        suffix = 2
        while new_name in used_names:
            new_name = f"{base_name}_{suffix}"
            suffix += 1
        used_names.add(new_name)
        chain = model.add_chain(source_chain)
        chain.name = new_name
        for residue in chain:
            residue.subchain = f"{new_name}_{residue.subchain or original_name}"
            residue.entity_id = f"{new_name}_{residue.entity_id or 'unknown'}"
        rows.append(
            {
                "source_id": source_id,
                "source_file": source_file,
                "original_chain": source_chain.name,
                "overlay_chain": new_name,
                "residues": len(chain),
                "atoms": sum(len(residue) for residue in chain),
                "heme_residues": sum(residue.name == "HEM" for residue in chain),
            }
        )
    return rows


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    manifest = pl.read_csv(input_dir / "alignment_manifest.csv").sort("pdb_id")
    if manifest.height != config["expected_source_count"]:
        raise ValueError("Unexpected number of aligned source structures")

    overlay = gemmi.Structure()
    overlay.name = "CYP3A4_public_overlay"
    model = overlay.add_model(gemmi.Model("1"))
    used_names: set[str] = set()
    reference_path = input_dir / "afdb_reference.cif"
    reference = gemmi.read_structure(str(reference_path))
    rows = add_structure(model, reference, "AFDB", "afdb_reference.cif", used_names)
    reference_atoms = reference[0].count_atom_sites()
    total_atoms = reference_atoms
    total_hemes = 0

    for source in manifest.iter_rows(named=True):
        pdb_id = source["pdb_id"]
        alias = f"aligned_{pdb_id}.cif"
        path = input_dir / alias
        if sha256(path) != source["aligned_sha256"]:
            raise ValueError(f"Aligned CIF hash mismatch for {pdb_id}")
        structure = gemmi.read_structure(str(path))
        atom_count = structure[0].count_atom_sites()
        if atom_count != source["all_atom_sites"]:
            raise ValueError(f"Atom count mismatch for {pdb_id}")
        new_rows = add_structure(model, structure, pdb_id.upper(), alias, used_names)
        hemes = sum(row["heme_residues"] for row in new_rows)
        if hemes != source["heme_residues"]:
            raise ValueError(f"HEME count mismatch for {pdb_id}")
        rows.extend(new_rows)
        total_atoms += atom_count
        total_hemes += hemes

    output_path = output_dir / "cyp3a4_public_overlay_122.cif"
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
    overlay.make_mmcif_document(groups).write_file(str(output_path))
    reloaded = gemmi.read_structure(str(output_path))
    if len(reloaded) != 1 or reloaded[0].count_atom_sites() != total_atoms:
        raise ValueError("Combined CIF atom count failed round-trip validation")
    reloaded_hemes = sum(
        residue.name == "HEM" for chain in reloaded[0] for residue in chain
    )
    if reloaded_hemes != total_hemes or len(reloaded[0]) != len(rows):
        raise ValueError("Combined CIF chain/HEME count failed round-trip validation")
    if len({chain.name for chain in reloaded[0]}) != len(rows):
        raise ValueError("Combined CIF has duplicate chain names")

    pl.DataFrame(rows).write_csv(output_dir / "chain_map.csv")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "alignment_run_id": config["alignment_run_id"],
        "reference_id": config["reference_id"],
        "source_structures": manifest.height,
        "coordinate_models": len(reloaded),
        "chains": len(reloaded[0]),
        "atoms": total_atoms,
        "heme_residues": total_hemes,
        "file": output_path.name,
        "bytes": output_path.stat().st_size,
        "sha256": sha256(output_path),
        "representation": "All source chains and AFDB reference in one coordinate model",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
