# Every public CYP3A4 chain independently superposed onto AFDB

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

After mmCIF read-back, every source protein was aligned again to AFDB as an
independent check. All 142 residual rigid transforms were identity to the
reported precision. The trimmed core RMSD spans 0.62–1.65 Å (median 0.83 Å).
This verifies protein-core placement, not that all ligand poses coincide or
that the ensemble represents a single biological complex.
