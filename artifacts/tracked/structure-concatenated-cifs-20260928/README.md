# Concatenated aligned CYP3A4 mmCIFs

Recorded run: `shared/20260928T083019116710Z-45e44833`.
The requested single file is `outputs/cyp3a4_aligned_122.cif` in that run,
73,788,955 bytes, SHA-256
`fee2c1ce0e58177e2eb019e8ec9aba23ed95399e3296868a372e09957f6813ef`.
The payload remains local; this directory tracks its summary and block map.

This file is a plain concatenation of the 122 previously aligned mmCIF files
in PDB ID order. Each original `data_` block, coordinate, chain ID, HEM, and
ligand remains unchanged. Only a newline separator is added between files.
Gemmi parsed the result as 122 data blocks with 531,617 atom sites in total,
matching the input blocks. `block_map.csv` records each PDB ID, block name,
source hash, byte length, and start offset.

The [earlier combined overlay](../structure-combined-overlay-20260928/README.md)
instead placed all atoms into one coordinate model and renamed chains; it is
retained as a separate visualization artifact, not this plain concatenation.
