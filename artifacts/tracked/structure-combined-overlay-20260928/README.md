# Combined CYP3A4 public structure overlay

This earlier overlay fitted only one CYP3A4 chain per PDB entry; 20 additional
protein chains retained their original offsets. For an overlay in which every
source chain is independently fitted to AFDB, use the
[corrected all-chain file](../structure-all-chains-afdb-20260928/README.md).

Recorded run: `shared/20260928T082350706757Z-7e279917`.
The single combined mmCIF is
`outputs/cyp3a4_public_overlay_122.cif` in that run. It is 87,571,855 bytes
with SHA-256
`9457f5535e784baa604318433bd7522194a6388aee678463a8e470bdf6b00d00`.
The coordinate payload remains local; this directory tracks its compact
summary and chain map.

The file places the unmodified AFDB P08684 reference and all 122 previously
transformed public source structures into **one coordinate model**. It contains
535,651 atoms in 143 uniquely named chains,
including 142 HEM residues. Multiple protein copies, original bound ligands,
and other deposited atoms from each source are retained. `chain_map.csv` maps
each overlay chain back to its PDB entry and original chain ID; the AFDB
reference uses source ID `AFDB`. The file was read back with Gemmi to verify
its model, chain, atom, and HEM counts.

All source coordinates were already aligned to the AFDB frame in
`shared/20260927T160227340394Z-eb38b998`. The combined file is a visualization
overlay, not a single biological complex or a prediction for any blinded
ligand. The superposed copies overlap extensively by design.
