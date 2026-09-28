# Public CYP3A4 structure source pool

Recorded collection: `shared/20260927T154642066309Z-7c506d38` (2026-09-28 JST).
This is a source pool for the 20 blinded CYP3A4 structure targets, not a curated
ligand-pose training set.

The run called TogoMCP 2.7.4 `search_pdb_entity` against PDBj with query
`CYP3A4`. It preserved all 202 returned hits in `outputs/togo_pdbj_search.json`
(SHA-256 `77f8ff54cf1ed6075317cfe61df44a5c74d55ce030e5127eb590006a4a09aa60`).
Free-text search included unrelated proteins, so the collector selected hits
whose reported ligands field contains both `Cytochrome P450 3A4` and
`PROTOPORPHYRIN IX CONTAINING FE`, case-insensitively. It downloaded 122
experimental mmCIF files from the RCSB file endpoint, totaling 102,836,041
bytes. All 122 selected hits are X-ray structures; resolutions span 1.7–3.4 Å
(median 2.55 Å). The mmCIF files remain in the recorded run because the pool is
too large for ordinary Git storage.

`source_manifest.csv` records each PDB ID, method, resolution, summary of
reported components, citation DOI, retrieval URL, local run-relative path, size,
and SHA-256. The run also snapshots the official 20-molecule CSV from dataset
revision `146282d2ce1532f1803679a834a9bec9dac1e482` and CYP3A4 FASTA from
tutorial revision `832ae6d`.

Selection by search-result text does not establish that a structure has a bound
drug-like ligand in the CYP3A4 active site, a comparable protein construct, or
a useful pose for a blinded target. Those checks and ligand similarity ranking
are separate curation work. Do not use the 122 files as labelled complexes
without that review.
