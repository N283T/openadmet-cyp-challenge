# CYP3A4 structures in the AFDB coordinate frame

AFDB snapshot run: `shared/20260927T160005116456Z-de79ae45`.
Alignment run: `shared/20260927T160227340394Z-eb38b998`.
The input source pool is
[`shared/20260927T154642066309Z-7c506d38`](../structure-cyp3a4-public-sources-20260928/README.md).

The reference is the [AlphaFold DB CYP3A4 model](https://alphafold.ebi.ac.uk/entry/P08684),
`AF-P08684-F1-model_v6`, for UniProt `P08684`. Its 503-residue sequence matches
the official tutorial FASTA. The run preserves the AFDB API response, original
PDB and mmCIF, URLs, and SHA-256 hashes. The AFDB model does not contain HEM.

Gemmi 0.7.5 aligned each public structure's protein sequence to the AFDB
protein and superposed Cα atoms with five outlier-trimming cycles. It selected
the source chain with the largest matched Cα count, then applied that chain's
rigid transform to **all atoms** in the source structure, including HEM and
bound ligands. Each transformed complex is in the alignment run's
`outputs/aligned/{pdb_id}.cif`; `outputs/afdb_reference.cif` is the unmodified
reference. Transformed mmCIFs omit crystal cell and symmetry metadata because
their coordinates are no longer in the original crystallographic frame. Files
with multiple protein copies retain the full source asymmetric unit; the
`aligned_chain` column identifies the copy placed onto the AFDB reference.

All 122 inputs aligned and retain their atom and HEM counts after mmCIF
round-trip validation. Matched Cα counts are at least 319 (median 343). The
outlier-trimmed core RMSD spans 0.63–1.65 Å, median 0.84 Å. `alignment_manifest.csv`
records the chain, matched count, RMSD, HEM count, input and output hashes, and
run-relative output path. As [Gemmi documents](https://gemmi.readthedocs.io/en/stable/analysis.html#superposition),
its RMSD convention divides by 3N; the reported `core_rmsd_angstrom` converts
this by multiplying by √3. It measures the selected structural core and does
not measure flexible-loop or ligand-pose agreement.

This is a coordinate-frame collection for inspection and future curation, not
a predicted structure for any of the 20 blinded ligands. Ligand identity and
active-site pose relevance remain to be assessed.
