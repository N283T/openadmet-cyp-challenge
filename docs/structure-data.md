# CYP structure and MSA source data

Four CYP targets are now prepared for prospective `boltz-trunk` models. The
[compact four-target record](../artifacts/tracked/cyp-structure-prep-20260928/README.md)
links immutable runs and contains source manifests, chain maps, ligand audits,
contact tables, and pocket definitions.

The initial trunk input and runtime settings are recorded separately in
[Boltz trunk input configuration](boltz-trunk.md). This route uses fixed HEM,
the extended pocket, MSA, and no template or affinity head.

| Target | Blocks | PDB entries | Retained HEM / ligands | Holo chains | Core / extended residues | MSA rows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CYP1A2 / P05177 | 2 | 1 | 1 / 1 | 1 | 20 / 20 | 12,515 |
| CYP2C9 / P11712 | 26 | 15 | 25 / 25 | 23 | 21 / 33 | 15,696 |
| CYP2D6 / P10635 | 49 | 14 | 48 / 54 | 44 | 18 / 33 | 16,017 |
| CYP3A4 / P08684 | 143 | 122 | 142 / 138 | 129 | 11 / 33 | 14,300 |

Each new ensemble is one CIF with multiple `data_` blocks for convenient
viewer selection. Every experimental block retains one AFDB-aligned CYP chain,
its corresponding HEM or older HEC/`HEME C` cofactor, and accepted active-site
co-crystal ligands. HOH, configured crystallization factors, and unrelated
non-polymers are removed. Accepted ligands have at least six heavy atoms,
contact protein within 4.5 Å, lie within 12 Å of a heme iron, and are not in
the configured artifact list. Residues are mapped by sequence alignment to
canonical 1-based UniProt numbering. The AFDB reference blocks have no HEM or
ligand because the source AFDB models contain neither.

The A3Ms are from pinned LiteFold revision
`632cfc87f567e429570897c3addd80665a697f3e`; they are not AFDB's own model-input
MSAs. Query sequences match their canonical AFDB model exactly. CYP1A2 has only
one exact CYP1A2/HEME structure in the current TogoMCP result, so its pocket is
not yet a frequency-stable structural consensus.

## Earlier CYP3A4 structure-track work

The blinded structure track has 20 ligand targets for human CYP3A4 (UniProt
P08684). The official tutorial FASTA and the AFDB v6 model use the same
503-residue sequence. The source data below support inspection and future
modeling; none is a predicted complex for a blinded ligand.

| Data | Recorded run | Compact record | Local payload |
| --- | --- | --- | --- |
| 122 public CYP3A4/HEM X-ray mmCIFs | `shared/20260927T154642066309Z-7c506d38` | [Source pool](../artifacts/tracked/structure-cyp3a4-public-sources-20260928/README.md) | `outputs/pdb/{pdb_id}.cif` |
| AFDB v6 full-length CYP3A4 model | `shared/20260927T160005116456Z-de79ae45` | [AFDB and alignment](../artifacts/tracked/structure-afdb-superposition-20260928/README.md) | `outputs/afdb_model.cif`, `outputs/afdb_model.pdb` |
| Public structures superposed onto AFDB | `shared/20260927T160227340394Z-eb38b998` | [Alignment manifest](../artifacts/tracked/structure-afdb-superposition-20260928/alignment_manifest.csv) | `outputs/aligned/{pdb_id}.cif` |
| **One CIF; AFDB and each of 142 fitted chains in separate data blocks** | `shared/20260928T090028397804Z-e19d26c5` | [Current multi-block CIF](../artifacts/tracked/structure-multiblock-superposed-20260928/README.md) | `outputs/cyp3a4_all_chains_143_data_blocks.cif` |
| Earlier all-chain fit in one data block | `shared/20260928T085009144228Z-c9f9fe2d` | [One-model source](../artifacts/tracked/structure-all-chains-afdb-20260928/README.md) | `outputs/cyp3a4_all_142_chains_on_afdb.cif` |
| AFDB plus all aligned CIFs, unchanged separate data blocks | `shared/20260928T084148280156Z-34becdd5` | [Complete concatenated CIF](../artifacts/tracked/structure-concatenated-afdb-cifs-20260928/README.md) | `outputs/cyp3a4_afdb_aligned_123.cif` |
| Earlier aligned CIFs concatenated without AFDB | `shared/20260928T083019116710Z-45e44833` | [122-block CIF](../artifacts/tracked/structure-concatenated-cifs-20260928/README.md) | `outputs/cyp3a4_aligned_122.cif` |
| Earlier one-model overlay with 20 extra chains not individually fitted | `shared/20260928T082350706757Z-7e279917` | [Earlier overlay](../artifacts/tracked/structure-combined-overlay-20260928/README.md) | `outputs/cyp3a4_public_overlay_122.cif` |
| Coordinate coverage against P08684 | `shared/20260927T161758979224Z-4d82a646` | [Coverage tables](../artifacts/tracked/structure-public-coverage-20260928/README.md) | `outputs/structure_coverage.csv`, `outputs/residue_coverage.csv` |
| Alternative full-length P08684 MSA | `shared/20260927T165502581117Z-7ca2b7c0` | [MSA provenance](../artifacts/tracked/structure-cyp3a4-msa-20260928/README.md) | `outputs/P08684.a3m` |

Local payload paths are relative to `artifacts/runs/shared/<run_id>/`.
External coordinates and the A3M remain in their recorded runs; compact
provenance and analysis tables are tracked in Git. The public structure search
is a source pool, not a curated set of ligand poses. The superposition places
protein, HEM, and ligands in a common coordinate frame, but does not establish
pose suitability for any target.
The current CIF has one file with 143 separate `data_` blocks. Each block
contains exactly one chain, and all 142 source chains were fitted independently
to the AFDB reference. This supports selecting blocks in a viewer without
loading one giant coordinate model.
The earlier all-chain overlay places the AFDB reference and all 142 source
chains in one coordinate model after fitting each source chain independently.
The complete concatenation instead preserves AFDB and each earlier aligned CIF
as separate `data_` blocks; many viewers show only one block at a time. The
earlier one-model overlay retained offsets for 20 extra protein chains.
None of these files is a single biological complex.

The P08684 A3M contains 14,300 aligned sequences across 503 query columns,
with a query matching the official full-length sequence. It was generated by
ColabFold/MMseqs2 and is **not** the MSA used for the AFDB v6 prediction.
AFDB advertises [its own v6 MSA](https://alphafold.ebi.ac.uk/files/msa/AF-P08684-F1-msa_v6.a3m),
but that URL returned HTTP 403 on 2026-09-28. The failed retrieval is recorded
as `shared/20260927T164601418477Z-3e8c6d36`.

Sources and conclusions are recorded in [RL0192–RL0196](research/2026-09.md#rl0192--public-cyp3a4-heme-structure-source-pool-collected).
No structure-track submission has been prepared or uploaded.
