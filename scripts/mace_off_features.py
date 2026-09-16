"""Extract ligand-only MACE-OFF23 atomic embeddings from deterministic conformers."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from ase import Atoms
from mace.calculators import MACECalculator
from rdkit import Chem
from rdkit.Chem import AllChem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
model_path = inp / cfg["model"]
model_digest = hashlib.sha256(model_path.read_bytes()).hexdigest()
if model_digest != cfg["model_sha256"]:
    raise ValueError(f"MACE model hash mismatch: {model_digest}")

compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
if cfg["max_compounds"] is not None:
    compounds = compounds.head(cfg["max_compounds"])

calculator = MACECalculator(
    model_paths=str(model_path),
    device=cfg["device"],
    default_dtype=cfg["default_dtype"],
)
supported = set(cfg["supported_atomic_numbers"])
embedding_rows = []
property_rows = []
failure_rows = []
embedding_width = None

for row_index, record in enumerate(compounds.iter_rows(named=True)):
    molecule = Chem.AddHs(Chem.MolFromSmiles(record["standardized_smiles"]))
    atomic_numbers = {atom.GetAtomicNum() for atom in molecule.GetAtoms()}
    if not atomic_numbers.issubset(supported):
        failure_rows.append(
            {
                "compound_id": record["compound_id"],
                "reason": "unsupported_elements",
                "detail": ",".join(map(str, sorted(atomic_numbers - supported))),
            }
        )
        continue
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = cfg["conformer_seed"]
    parameters.useRandomCoords = True
    status = AllChem.EmbedMolecule(molecule, parameters)
    if status != 0:
        failure_rows.append(
            {
                "compound_id": record["compound_id"],
                "reason": "conformer_embedding",
                "detail": str(status),
            }
        )
        continue
    conformer = molecule.GetConformer()
    ase_atoms = Atoms(
        numbers=[atom.GetAtomicNum() for atom in molecule.GetAtoms()],
        positions=np.asarray(conformer.GetPositions(), dtype=float),
        pbc=False,
    )
    try:
        descriptors = calculator.get_descriptors(ase_atoms, invariants_only=True)
        ase_atoms.calc = calculator
        energy = float(ase_atoms.get_potential_energy())
    except Exception as error:  # Preserve per-compound model failures without dropping rows.
        failure_rows.append(
            {
                "compound_id": record["compound_id"],
                "reason": type(error).__name__,
                "detail": str(error),
            }
        )
        continue
    if embedding_width is None:
        embedding_width = descriptors.shape[1]
    if descriptors.shape[1] != embedding_width or not np.isfinite(descriptors).all():
        raise ValueError("MACE descriptor shape or finiteness changed within the run")
    pooled = np.concatenate([descriptors.mean(axis=0), descriptors.max(axis=0)])
    embedding_rows.append(
        {
            "compound_id": record["compound_id"],
            **{f"mace_off_emb_{index}": float(value) for index, value in enumerate(pooled)},
        }
    )
    property_rows.append(
        {
            "compound_id": record["compound_id"],
            "mace_off_energy": energy,
            "mace_off_energy_per_atom": energy / molecule.GetNumAtoms(),
            "mace_off_atom_count_with_h": molecule.GetNumAtoms(),
        }
    )
    if (row_index + 1) % 100 == 0:
        print(f"MACE-OFF attempted {row_index + 1}/{compounds.height}", flush=True)

failures = pl.DataFrame(failure_rows) if failure_rows else pl.DataFrame(
    schema={"compound_id": pl.String, "reason": pl.String, "detail": pl.String}
)
failures.write_parquet(out / "mace-off-failures.parquet")
embeddings = pl.DataFrame(embedding_rows)
properties = pl.DataFrame(property_rows)
embeddings.write_parquet(out / "mace-off-embeddings.parquet")
properties.write_parquet(out / "mace-off-properties.parquet")

summary = {
    "code_commit": cfg["code_commit"],
    "model_commit": cfg["model_commit"],
    "model_sha256": model_digest,
    "model": "MACE-OFF23 small",
    "attempted_compounds": compounds.height,
    "successful_compounds": embeddings.height,
    "failed_compounds": failures.height,
    "atomic_embedding_width": embedding_width,
    "pooled_embedding_width": 0 if embedding_width is None else 2 * embedding_width,
    "pooling": "mean and componentwise maximum over explicit-H atom descriptors",
    "conformer": "RDKit ETKDGv3, one unoptimized conformer",
    "limitations": [
        "MACE-OFF23 is trained on neutral closed-shell organic chemistry.",
        "The embedding is ligand-only and is not a CYP reaction or TDI model.",
        "Unsupported elements and conformer failures remain explicitly missing.",
    ],
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
