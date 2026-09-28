# AFDB reference and aligned CYP3A4 CIFs in one file

Recorded run: `shared/20260928T084148280156Z-34becdd5`.
The current complete file is `outputs/cyp3a4_afdb_aligned_123.cif` in that run,
74,262,970 bytes, SHA-256
`fe36c877e1b74cdbbd01c52130bbc49685e9e939357253d8768653fb0222c4ab`.
The coordinate payload remains local; this directory tracks the summary and
block map.

The first unchanged `data_` block is the full-length AFDB reference
`AF-P08684-F1` (4,034 atoms). The next 122 blocks are the previously aligned
public experimental CYP3A4 structures, each with its original chain IDs,
coordinates, HEM, and ligands unchanged. Only newline separators were added.
Gemmi parsed all 123 blocks and verified 535,651 atom sites in total. The
`block_map.csv` lists each source, block name, original SHA-256, byte length,
and starting byte offset.

The [earlier 122-block file](../structure-concatenated-cifs-20260928/README.md)
omitted AFDB and is retained as an earlier result. The separate
[one-model overlay](../structure-combined-overlay-20260928/README.md) copied
atoms into renamed chains and serves a different visualization purpose.
