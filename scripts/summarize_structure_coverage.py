"""Map public CYP3A4 coordinate coverage to UniProt P08684 positions."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import polars as pl


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    reference = gemmi.read_structure(str(input_dir / "afdb_model.cif"))[0][0].get_polymer()
    reference_residues = reference.extract_sequence()
    reference_sequence = gemmi.one_letter_code(reference_residues)
    if len(reference_sequence) != config["reference_length"]:
        raise ValueError("Unexpected AFDB reference length")

    sources = pl.read_csv(input_dir / "source_manifest.csv")
    alignments = pl.read_csv(input_dir / "alignment_manifest.csv")
    if sources["pdb_id"].to_list() != alignments["pdb_id"].to_list():
        raise ValueError("Source and alignment manifests have different PDB IDs")
    counts = [0] * len(reference_sequence)
    rows = []
    for source_row, alignment_row in zip(
        sources.iter_rows(named=True), alignments.iter_rows(named=True), strict=True
    ):
        pdb_id = source_row["pdb_id"]
        path = input_dir / f"pdb_{pdb_id}.cif"
        if hashlib.sha256(path.read_bytes()).hexdigest() != source_row["sha256"]:
            raise ValueError(f"Source hash mismatch for {pdb_id}")
        polymer = gemmi.read_structure(str(path))[0][alignment_row["aligned_chain"]].get_polymer()
        alignment = gemmi.align_sequence_to_polymer(
            reference_residues,
            polymer,
            gemmi.PolymerType.PeptideL,
            gemmi.AlignmentScoring("b"),
        )
        gapped_reference = alignment.add_gaps(reference_sequence, 1)
        gapped_source = alignment.add_gaps(gemmi.one_letter_code(polymer.extract_sequence()), 2)
        position = 0
        covered = []
        for reference_letter, source_letter in zip(
            gapped_reference, gapped_source, strict=True
        ):
            if reference_letter != "-":
                position += 1
                if source_letter != "-":
                    counts[position - 1] += 1
                    covered.append(position)
        if position != len(reference_sequence) or not covered:
            raise ValueError(f"Invalid sequence mapping for {pdb_id}")
        first, last = covered[0], covered[-1]
        rows.append(
            {
                "pdb_id": pdb_id,
                "aligned_chain": alignment_row["aligned_chain"],
                "coordinate_polymer_residues": len(polymer),
                "mapped_reference_residues": len(covered),
                "reference_coverage_fraction": len(covered) / len(reference_sequence),
                "sequence_identity_fraction": alignment.match_count / len(covered),
                "first_mapped_position": first,
                "last_mapped_position": last,
                "missing_n_terminus": first - 1,
                "missing_c_terminus": len(reference_sequence) - last,
                "missing_between_termini": (last - first + 1) - len(covered),
            }
        )
    structures = pl.DataFrame(rows).sort("pdb_id")
    structures.write_csv(output_dir / "structure_coverage.csv")
    residues = pl.DataFrame(
        {
            "uniprot_position": list(range(1, len(reference_sequence) + 1)),
            "amino_acid": list(reference_sequence),
            "structures_with_coordinates": counts,
        }
    ).with_columns((pl.col("structures_with_coordinates") / sources.height).alias("fraction"))
    residues.write_csv(output_dir / "residue_coverage.csv")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "source_run_id": config["source_run_id"],
        "afdb_run_id": config["afdb_run_id"],
        "alignment_run_id": config["alignment_run_id"],
        "structure_count": structures.height,
        "reference_residues": len(reference_sequence),
        "mapped_min": structures["mapped_reference_residues"].min(),
        "mapped_median": structures["mapped_reference_residues"].median(),
        "mapped_max": structures["mapped_reference_residues"].max(),
        "first_mapped_median": structures["first_mapped_position"].median(),
        "last_mapped_median": structures["last_mapped_position"].median(),
        "structures_at_least_90_percent": structures.filter(
            pl.col("reference_coverage_fraction") >= 0.9
        ).height,
        "positions_in_all_structures": residues.filter(
            pl.col("structures_with_coordinates") == structures.height
        ).height,
        "method": "Gemmi sequence-to-polymer alignment with BLOSUM62; coordinate residues only",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
