"""Independently superpose every public CYP3A4 chain onto one AFDB model."""

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


def count_hemes(model: gemmi.Model) -> int:
    return sum(residue.name == "HEM" for chain in model for residue in chain)


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
    sources = pl.read_csv(input_dir / "source_manifest.csv").sort("pdb_id")
    if sources.height != config["expected_source_count"]:
        raise ValueError("Unexpected number of source structures")

    reference = gemmi.read_structure(str(input_dir / "afdb_model.cif"))
    if len(reference) != 1 or len(reference[0]) != 1:
        raise ValueError("Expected one AFDB reference chain")
    reference_polymer = reference[0][0].get_polymer()
    if len(reference_polymer) != 503:
        raise ValueError("Expected full-length 503-residue AFDB CYP3A4")
    protein_type = reference_polymer.check_polymer_type()

    overlay = gemmi.Structure()
    overlay.name = "CYP3A4_all_chains_on_AFDB"
    model = overlay.add_model(gemmi.Model("1"))
    afdb_chain = model.add_chain(reference[0][0])
    afdb_chain.name = "AFDB_A"
    for residue in afdb_chain:
        residue.subchain = f"AFDB_{residue.subchain or 'A'}"
        residue.entity_id = f"AFDB_{residue.entity_id or '1'}"
    total_atoms = reference[0].count_atom_sites()
    total_hemes = 0
    rows = []
    for source in sources.iter_rows(named=True):
        pdb_id = source["pdb_id"]
        source_path = input_dir / f"pdb_{pdb_id}.cif"
        if sha256(source_path) != source["sha256"]:
            raise ValueError(f"Source hash mismatch for {pdb_id}")
        structure = gemmi.read_structure(str(source_path))
        if len(structure) != 1:
            raise ValueError(f"Expected one model in {pdb_id}")
        total_atoms += structure[0].count_atom_sites()
        total_hemes += count_hemes(structure[0])
        for chain in structure[0]:
            polymer = chain.get_polymer()
            if len(polymer) < config["min_matched_ca"]:
                raise ValueError(f"Unrecognized non-CYP3A4 chain in {pdb_id}:{chain.name}")
            sup = gemmi.calculate_superposition(
                reference_polymer,
                polymer,
                protein_type,
                gemmi.SupSelect.CaP,
                trim_cycles=config["trim_cycles"],
                trim_cutoff=config["trim_cutoff"],
            )
            core_rmsd = sup.rmsd * math.sqrt(3)
            if sup.count < config["min_matched_ca"] or core_rmsd > config["max_core_rmsd_angstrom"]:
                raise ValueError(f"Poor CYP3A4 fit for {pdb_id}:{chain.name}")

            single = gemmi.Model("1")
            moved = single.add_chain(chain)
            single.transform_pos_and_adp(sup.transform)
            new_name = f"{pdb_id.upper()}_{chain.name or 'blank'}"
            if model.find_chain(new_name) is not None:
                raise ValueError(f"Duplicate combined chain name: {new_name}")
            added = model.add_chain(moved)
            added.name = new_name
            for residue in added:
                residue.subchain = f"{new_name}_{residue.subchain or chain.name}"
                residue.entity_id = f"{new_name}_{residue.entity_id or 'unknown'}"
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "original_chain": chain.name,
                    "overlay_chain": new_name,
                    "matched_ca": sup.count,
                    "core_rmsd_angstrom": core_rmsd,
                    "polymer_residues": len(polymer),
                    "heme_residues": sum(residue.name == "HEM" for residue in added),
                    "atoms": sum(len(residue) for residue in added),
                    "source_sha256": source["sha256"],
                }
            )

    if len(rows) != config["expected_source_chains"]:
        raise ValueError(f"Expected 142 source chains, got {len(rows)}")
    output_path = output_dir / "cyp3a4_all_142_chains_on_afdb.cif"
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
    checked = gemmi.read_structure(str(output_path))
    if len(checked) != 1 or len(checked[0]) != len(rows) + 1:
        raise ValueError("Combined mmCIF model/chain count failed round-trip validation")
    if checked[0].count_atom_sites() != total_atoms or count_hemes(checked[0]) != total_hemes:
        raise ValueError("Combined mmCIF atom/HEME count failed round-trip validation")

    for row in rows:
        chain = checked[0][row["overlay_chain"]]
        fit = gemmi.calculate_superposition(
            reference_polymer,
            chain.get_polymer(),
            protein_type,
            gemmi.SupSelect.CaP,
            trim_cycles=config["trim_cycles"],
            trim_cutoff=config["trim_cutoff"],
        )
        translation, rotation = residual_transform(fit.transform)
        if (
            translation > config["max_residual_translation_angstrom"]
            or rotation > config["max_residual_rotation_element"]
        ):
            raise ValueError(f"Chain is not superposed onto AFDB: {row['overlay_chain']}")
        row["residual_translation_angstrom"] = translation
        row["residual_rotation_element"] = rotation

    report = pl.DataFrame(rows).sort("pdb_id", "original_chain")
    report.write_csv(output_dir / "chain_alignment.csv")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "source_run_id": config["source_run_id"],
        "afdb_run_id": config["afdb_run_id"],
        "source_structures": sources.height,
        "independently_superposed_source_chains": report.height,
        "coordinate_models": len(checked),
        "chains_including_afdb": len(checked[0]),
        "atoms": total_atoms,
        "heme_residues": total_hemes,
        "core_rmsd_min_angstrom": report["core_rmsd_angstrom"].min(),
        "core_rmsd_median_angstrom": report["core_rmsd_angstrom"].median(),
        "core_rmsd_max_angstrom": report["core_rmsd_angstrom"].max(),
        "max_residual_translation_angstrom": report["residual_translation_angstrom"].max(),
        "max_residual_rotation_element": report["residual_rotation_element"].max(),
        "file": output_path.name,
        "bytes": output_path.stat().st_size,
        "sha256": sha256(output_path),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
