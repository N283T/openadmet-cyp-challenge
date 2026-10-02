# CYP structure and MSA preparation

Compact, Git-tracked records for the four CYP targets prepared for prospective
`boltz-trunk` models. Large source mmCIFs, AFDB coordinates, A3Ms, and combined
ensemble CIFs remain in their immutable recorded runs.

| Target | Accession | Reference + MSA run | TogoMCP/PDB run | Ensemble + pocket run |
| --- | --- | --- | --- | --- |
| CYP1A2 | P05177 | `shared/20260928T111653134056Z-84ffca6a` | `shared/20260928T111808486836Z-1dc3f7ab` | `shared/20260928T113128124677Z-66cecb8c` |
| CYP2C9 | P11712 | `shared/20260928T111758347957Z-0357cb7c` | `shared/20260928T112422815798Z-83412439` | `shared/20260928T113128145767Z-366f1503` |
| CYP2D6 | P10635 | `shared/20260928T111758358079Z-0d42000d` | `shared/20260928T111808525822Z-5b7c2291` | `shared/20260928T113128150604Z-acf3b0b7` |
| CYP3A4 | P08684 | `shared/20260928T111713518371Z-c156bd1a` | `shared/20260928T111808507687Z-e19a6977` | `shared/20260928T113128171818Z-ca17bfa3` |

Each target directory contains the reference, source, and ensemble summaries;
the TogoMCP-derived source manifest; the ensemble block map; the complete
ligand inclusion/exclusion audit; individual ligand contacts; and pocket
definitions in CSV and JSON.

The combined CIF for each target has one `data_` block for the full-length AFDB
reference followed by one block per accepted experimental chain. Each
experimental block retains exactly one selected CYP protein chain, its
corresponding HEM/HEC cofactor, and every accepted active-site co-crystal
ligand. HOH, configured crystallization artifacts, and unrelated non-polymers
are removed. The AFDB block has no HEM or ligand because neither is present in
the source AFDB model.

Pocket residue numbering is the 1-based canonical UniProt/AFDB sequence. A
candidate pocket ligand must not be water, HEM/HEC, or a configured crystallization
artifact; must have at least six heavy atoms; must contact the protein within
4.5 Å; and must lie within 12 Å of a HEM iron. `core_residues` contact at least
50% of accepted holo chains. `extended_residues` include the core plus residues
contacting at least 10% and at least two holo chains. CYP1A2 has only one public
structure in this exact TogoMCP selection, so its 20-residue pocket is a
single-structure definition and must not be interpreted as a frequency-stable
consensus.

The LiteFold A3Ms come from pinned dataset revision
`632cfc87f567e429570897c3addd80665a697f3e`. Their queries and the AFDB
canonical models were independently required to match the same full-length
sequence. They are alternative ColabFold/MMseqs2 MSAs, not AFDB's private
prediction input MSAs.
