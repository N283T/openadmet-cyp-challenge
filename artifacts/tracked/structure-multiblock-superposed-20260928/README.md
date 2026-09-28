# All fitted CYP3A4 chains in one multi-block CIF

Recorded run: `shared/20260928T090028397804Z-e19d26c5`.
The requested file is `outputs/cyp3a4_all_chains_143_data_blocks.cif` in that
run. It is **one CIF file** containing 143 separate `data_` blocks: the AFDB
reference first, then one independently fitted CYP3A4 chain per block.
No individual chain CIFs or ZIP are produced. The 72,281,870-byte coordinate
file remains local; its SHA-256 is
`455ee7c1f2665af68bdf3f8b9f85ad118ee2d0d637e808329d0f4a1f0b0fbf73`.

Every source block retains its original chain ID within that block, with its
HEM, ligands, and other chain atoms. `block_map.csv` records the PDB ID,
original chain, block name, atom count, and alignment checks. Gemmi read back
the complete file as 143 blocks, each with exactly one chain: 535,651 atoms
and 142 HEM residues in total. Reconstructing polymer entities after parsing
and aligning each block to AFDB gave residual translations below
4 × 10⁻⁹ Å and rotation-matrix element differences below 3 × 10⁻¹⁰.

The [previous one-model file](../structure-all-chains-afdb-20260928/README.md)
puts all 143 chains in one data block, which is unsuitable for the requested
viewer workflow. This multi-block file is the current deliverable.
