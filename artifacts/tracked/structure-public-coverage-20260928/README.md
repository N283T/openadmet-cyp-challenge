# Coordinate coverage of the public CYP3A4 structures

Recorded run: `shared/20260927T161758979224Z-4d82a646`.
The run maps the selected CYP3A4 protein chain in each of the 122 public
structures to the 503-residue UniProt P08684 sequence represented by the AFDB
model. Gemmi aligns the full reference sequence to coordinate-bearing polymer
residues with BLOSUM62. `structure_coverage.csv` gives coverage for each PDB
entry; `residue_coverage.csv` gives the count for each UniProt position.

The number of mapped residues with coordinates is 415–470 per structure
(median 454), or 82.5–93.4% of the 503-residue reference (median 90.3%).
Seventy-four of 122 entries reach at least 90%; none has coordinates for all
503 residues. The first mapped position has median 28 (range 22–30), and the
last has median 496 (range 495–499). Positions 1–20 and 500–503 have no
coordinates in this collection. Internal variability is pronounced at 263–268
and 281–286; positions 282–285 appear in only one structure each. In total,
387 reference positions have coordinates in every structure.

Coverage describes observed coordinates in the deposited files, not the
sequence length of the expressed construct. The 319–425 Cα atoms used for
the separately reported trimmed superposition are a smaller structural core,
not the full coordinate coverage. These entries are public X-ray structures;
the blinded 20 cryo-EM target structures remain unavailable.
