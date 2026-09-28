# Every public CYP3A4 chain independently superposed onto AFDB

This earlier file has all chains in one data block. For the requested **one
CIF with one chain per `data_` block**, use the
[143-block file](../structure-multiblock-superposed-20260928/README.md).

Recorded run: `shared/20260928T085009144228Z-c9f9fe2d`.
The corrected simultaneous-display file is
`outputs/cyp3a4_all_142_chains_on_afdb.cif` in that run. It is 87,549,404
bytes, SHA-256
`074d7b221b6d41fdfea2bad41c10c7c6bbeac3932e543dca4ea75b6724f0da82`.
The large coordinate payload remains local; this directory tracks its summary
and per-chain alignment table.

All 142 CYP3A4 chains from the 122 public experimental entries were **each**
aligned to the full-length AFDB P08684 reference using their own Cα
superposition. Every atom belonging to a source chain, including its HEM,
ligands, and solvent, received that chain's transform. The combined mmCIF has
one coordinate model with 143 chains (142 sources plus AFDB), 535,651 atoms,
and 142 HEM residues. Unique PDB-prefixed chain names make simultaneous
display possible; `chain_alignment.csv` maps them to original PDB IDs and
chain IDs.

The original run's residual-transform check did not reconstruct polymer
entities after mmCIF parsing, so its reported zero residuals were not a valid
independent check. The subsequent 143-block run reconstructed those entities
and verified all 142 fits. The trimmed core RMSD spans 0.62–1.65 Å
(median 0.83 Å).
This verifies protein-core placement, not that all ligand poses coincide or
that the ensemble represents a single biological complex.
