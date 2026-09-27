"""Superpose full CYP3A4/HEM structures onto the AFDB P08684 model."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import polars as pl

MMCIF_GROUPS = gemmi.MmcifOutputGroups(
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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def heme_count(structure: gemmi.Structure) -> int:
    return sum(residue.name == "HEM" for chain in structure[0] for residue in chain)


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    manifest = pl.read_csv(input_dir / "source_manifest.csv")
    reference_path = input_dir / "afdb_model.cif"
    reference = gemmi.read_structure(str(reference_path))
    if len(reference) != 1 or len(reference[0]) != 1:
        raise ValueError("Expected a monomeric AFDB CYP3A4 reference")
    reference_polymer = reference[0][0].get_polymer()
    protein_type = reference_polymer.check_polymer_type()
    if len(reference_polymer) != 503:
        raise ValueError("AFDB reference is not the expected 503-residue CYP3A4 model")
    shutil.copyfile(reference_path, output_dir / "afdb_reference.cif")

    aligned_dir = output_dir / "aligned"
    aligned_dir.mkdir()
    rows = []
    for source_row in manifest.iter_rows(named=True):
        pdb_id = source_row["pdb_id"]
        source_path = input_dir / f"pdb_{pdb_id}.cif"
        if sha256(source_path) != source_row["sha256"]:
            raise ValueError(f"Source hash mismatch for {pdb_id}")
        structure = gemmi.read_structure(str(source_path))
        if not structure or heme_count(structure) == 0:
            raise ValueError(f"No HEM coordinates in {pdb_id}")
        atom_count = structure[0].count_atom_sites()
        candidates = []
        for chain in structure[0]:
            polymer = chain.get_polymer()
            if len(polymer) < config["min_matched_ca"]:
                continue
            try:
                sup = gemmi.calculate_superposition(
                    reference_polymer,
                    polymer,
                    protein_type,
                    gemmi.SupSelect.CaP,
                    trim_cycles=config["trim_cycles"],
                    trim_cutoff=config["trim_cutoff"],
                )
            except RuntimeError:
                continue
            if sup.count >= config["min_matched_ca"]:
                candidates.append((sup.count, sup.rmsd, chain.name, sup, len(polymer)))
        if not candidates:
            raise ValueError(f"No alignable CYP3A4 chain in {pdb_id}")
        count, gemmi_rmsd, chain_name, sup, polymer_length = min(
            candidates, key=lambda item: (-item[0], item[1], item[2])
        )
        rmsd = gemmi_rmsd * math.sqrt(3)
        if rmsd > config["max_core_rmsd_angstrom"]:
            raise ValueError(f"Excessive core RMSD for {pdb_id}: {rmsd:.2f} Å")
        for model in structure:
            model.transform_pos_and_adp(sup.transform)
        output_path = aligned_dir / f"{pdb_id}.cif"
        structure.make_mmcif_document(MMCIF_GROUPS).write_file(str(output_path))
        verified = gemmi.read_structure(str(output_path))
        if verified[0].count_atom_sites() != atom_count or heme_count(verified) == 0:
            raise ValueError(f"Atom/HEM round-trip validation failed for {pdb_id}")
        rows.append(
            {
                "pdb_id": pdb_id,
                "reference_id": "AF-P08684-F1-model_v6",
                "aligned_chain": chain_name,
                "chain_polymer_residues": polymer_length,
                "matched_ca": count,
                "core_rmsd_angstrom": rmsd,
                "gemmi_rmsd_angstrom": gemmi_rmsd,
                "heme_residues": heme_count(verified),
                "all_atom_sites": atom_count,
                "source_sha256": source_row["sha256"],
                "aligned_file": f"aligned/{pdb_id}.cif",
                "aligned_sha256": sha256(output_path),
            }
        )
        if len(rows) % 20 == 0 or len(rows) == manifest.height:
            print(f"Aligned {len(rows)}/{manifest.height} structures", flush=True)

    result = pl.DataFrame(rows).sort("pdb_id")
    result.write_csv(output_dir / "alignment_manifest.csv")
    result.write_parquet(output_dir / "alignment_manifest.parquet")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "source_run_id": config["source_run_id"],
        "afdb_run_id": config["afdb_run_id"],
        "reference_sha256": sha256(reference_path),
        "source_count": manifest.height,
        "aligned_count": result.height,
        "matched_ca_min": result["matched_ca"].min(),
        "matched_ca_median": result["matched_ca"].median(),
        "core_rmsd_min_angstrom": result["core_rmsd_angstrom"].min(),
        "core_rmsd_median_angstrom": result["core_rmsd_angstrom"].median(),
        "core_rmsd_max_angstrom": result["core_rmsd_angstrom"].max(),
        "method": "Gemmi sequence alignment + CA superposition, five outlier-trim cycles",
        "rmsd_note": "Reported core RMSD is Gemmi RMSD multiplied by sqrt(3)",
        "coordinate_note": "Aligned mmCIF omits crystallographic cell/symmetry metadata",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
