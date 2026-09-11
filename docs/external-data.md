# Source-preserving external CYP auxiliaries

The first external dataset contains **24,075 standardized compounds and 93,119
compound/readout labels** after removing every connectivity block present in the
challenge, including test and auxiliary-only molecules. It is an auxiliary dataset,
not extra rows of challenge pIC50. Classification optimization remains deferred.

## Original sources

### ChEMBL 36

Reuse the prior PXR project's local PostgreSQL mirror (`chembl_36`, socket `/tmp`,
port 5433), rather than downloading ChEMBL again. Its version table reports
ChEMBL_36, creation date 2025-07-28. The community reference used ChEMBL 37;
our acquisition and counts are not a reproduction of that version.

Raw export `20260911T052136820510Z-d4df2307` contains 155,013 IC50/AC50/Potency/Ki
activity rows across 33,018 ChEMBL molecules for human CYP1A2, CYP2C9, CYP2D6,
CYP3A4 and CYP2C19. Activity, assay, molecule, target and document identifiers,
relations, units, assay descriptions, quality flags, DOI/PubMed and source IDs
are retained. SQL is in `scripts/extract_chembl.sql`. The source transaction is
repeatable-read and read-only. CYP2C19 is auxiliary only, not a challenge endpoint.

### PubChem AID 1851

Original [AID 1851](https://pubchem.ncbi.nlm.nih.gov/bioassay/1851) description and
results are retrieved through [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest).
Run `20260911T052330309388Z-80507a6e` preserves 17,143 SIDs, 16,560 current CIDs,
85,715 SID/isoform observations, all original result fields, and CID-keyed isomeric
SMILES/InChIKeys. All 300-SID and 300-CID batches are audited for exact identifier
coverage. URLs, response hashes, status, sizes and retrieval times are retained.
The assay description reports AID version 1, revision 2.

Use the **deposited SMILES** attached to each assay SID for modeling; the independently
retrieved current CID structures are retained as an additional source snapshot.
Do not join observations by row index. Duplicate SIDs mapping to a standardized
compound are aggregated only after retaining their original records.

The original assay description explicitly identifies **P450-Glo/pro-luciferin
luminescence**, using recombinant human enzymes and a luciferase detection reagent.
This differs from the community code's characterization of the Veith panel as
fluorescent. The original description also says that inhibitors and some substrates
can reduce the signal. The auxiliary target is therefore the assay's continuous
response, not a direct-inhibition potency label inferred from it.

`Max_Response` is the maximum observed activity, usually at the highest tested
concentration. Every panel observation in this snapshot has a finite numeric
response, including inactive observations. No fabricated potency is assigned to
inactives. PubChem QC fields are retained; the snapshot's entries name supplier/QC
sources, rather than exposing a simple pass/fail criterion.

## Curation and decontamination

Adopted run: `20260911T052745709989Z-8792f422`; configuration:
`configs/external/curation-v2.json`; implementation: `scripts/curate_external.py`.

1. ChEMBL exact IC50 only, human single-protein target, confidence >=8, and assay
   types A/B/F (ADME/binding/functional). Require positive nM values, finite pChEMBL,
   no data-validity comment and no potential-duplicate flag. Verify pChEMBL agrees
   with `9 - log10(IC50_nM)` within rounding tolerance.
2. Exclude all ChEMBL PubChem imports (`assay_src_id=7`) from this first ChEMBL
   target family. The export confirms the same panel under `1851_1` through
   `1851_5`, each with 17,143 rows. Other imported panels are retained in the raw
   snapshot but also excluded from this initial comparison.
3. Conservatively exclude descriptions explicitly mentioning TDI, time-dependent,
   mechanism-based, irreversible or preincubation conditions. This deliberately
   loses some ordinary assay preincubations; it is not a claim that all
   preincubation means time-dependent inhibition. Refinement requires assay review.
4. Preserve censored relations and rejected measurements in the raw/selection
   snapshots; they do not become exact labels in this experiment. Ki is retained
   in the export but not mixed with IC50.
5. Standardize using the project's ChEMBL parent pipeline, retaining stereochemistry.
   Quarantine failed/excluded structures. Verify idempotent standardization and
   compound IDs for all challenge molecules.
6. Exclude external observations whose standardized InChIKey connectivity block
   matches **any** challenge molecule. This is stricter than fold-only filtering;
   challenge-overlap data is preserved, not used for auxiliary fitting.
7. ChEMBL: median per assay/compound, then median across assays per compound/isoform.
   PubChem: clip individual numeric responses to [-150,50], then median per
   standardized compound/isoform. Clipping is a declared modeling heuristic from
   the community approach, not a verified way to detect assay artifacts.

| Auxiliary isoform | ChEMBL IC50 labels | PubChem response labels |
|---|---:|---:|
| CYP1A2 | 1,234 | 16,244 |
| CYP2C9 | 2,440 | 16,244 |
| CYP2D6 | 2,362 | 16,244 |
| CYP3A4 | 4,588 | 16,244 |
| CYP2C19 | 1,275 | 16,244 |

These are distinct standardized compound/readout labels within a source, not
independent compounds summed across columns. Sixty-eight compounds occur in both
external sources. Their source-specific labels remain distinct targets.

The selected ChEMBL data spans pIC50 >=4; it does not directly fill the weak-potency
tail. PubChem contributes continuous responses for weak/inactive compounds without
pretending that the response is a measured pIC50.

The raw PubChem CYP2D6 panel contains 11,139 inactive observations, all with numeric
Max_Response. Clipping is not equally rare across isoforms: for CYP2C19 it affects
504 active, 51 inactive and 359 inconclusive observations in the raw panel. The
full `pubchem-response-coverage.csv` is retained; do not assume the community
article's clipping-rate claim holds for this snapshot.

Curation excludes 236 selected ChEMBL rows and 1,555 PubChem rows for challenge
connectivity overlap, and 3/330 rows for structure failures, respectively. Counts
are observation rows and include isoform repetitions. The adopted target and
compound tables are unique on their explicit keys.

## PostgreSQL lineage

Migration `012_external_auxiliaries.sql` adds external dataset, compound, target,
and source-link tables while sharing the existing molecular identities. Registration
run `20260911T053359215158Z-3e86b57c` atomically registers 24,075 compounds, 93,119
targets and 96,714 underlying source observations. Registration verifies every
input hash against the completed curation run. No prior source data is overwritten.
Raw descriptions, rejected rows, HTTP responses and full structure maps remain in
the immutable run artifacts.

## First model comparison

The five-fold comparison is complete; see the
[results and decision](analysis/2026-09-11-external-comparison.md). External training
improves the matched neural point estimate but does not beat the adopted TabPFN
baseline. The external dataset is retained; the neural recipe is not adopted.

Configuration: `configs/external/training-v1.json`. Use CheMeleon initialization,
regression only, and fixed UMAP outer folds. Both arms have the same 18-head
architecture/initial weights: four scored pIC50 heads, four challenge log2fc heads,
five ChEMBL IC50 heads and five PubChem response heads. External heads receive no
loss in the control arm.

- Challenge objective: mean standardized MSE across observed direct heads, plus
  0.3 times mean standardized MSE across log2fc heads.
- External arm adds 0.3 times the average of the two external source objectives.
  Each source objective averages observed heads. Each challenge batch is paired
  with 32 sampled ChEMBL and 32 sampled PubChem compounds; overlapping external
  sources retain both target families. This is a sampled schedule, not a full
  external-data epoch for each challenge epoch.
- Control: identical challenge batches, target scales, learning rate and stopping
  policy, with external weight zero. Adding external samples also changes the
  BatchNorm batch distribution and dropout RNG consumption; the comparison tests
  the overall auxiliary-training recipe, not a perfectly isolated gradient term.
- The outer assignments stay fixed. Connectivity overlaps with the outer fold
  are removed from all fitting roles; inner validation is split by connectivity
  within the remaining molecules. Both arms use identical memberships.
- Means/SDs for challenge labels are fitted only on inner-training molecules.
  External scales use the disjoint external dataset. Early stopping sees only
  inner-validation direct-regression MSE, never outer labels or LB scores.
- Adam at 1e-4; batch 64; seed 42 with fold-specific model seed; max 60 epochs,
  patience 10; FFN 256, dropout 0.1. Selected and final checkpoints, membership,
  history, keyed OOF and test predictions are retained for each run.

This is a matched direct-neural comparison, not a change to the adopted TabPFN 2.6
submission pipeline. The saved TabPFN OOF is included as a performance reference.
No calibration or submission is coupled to this experiment.

## Retained failures and revisions

- `20260911T052136820583Z-c16f503e`: original PubChem whole-assay GET exceeds the
  10,000-SID limit; the fetch helper also failed trying to recreate an immutable
  request manifest. Preserve the returned description/error, then use paginated
  retrieval and separate per-request metadata.
- `20260911T052556943962Z-f9e57655`: completed curation v1, superseded. Its B/F-only
  restriction omitted ADME assays. Audit `20260911T052656297028Z-020f1f5e` exposed
  the attrition, leading to v2. No model was trained on v1.
- `20260911T053107722245Z-7f8addde`: initial control run failed before training;
  a suffix glob confused `targets.parquet` with `external-targets.parquet`.
  Exact snapshot-basename resolution and a regression test fix this ambiguity.
- `20260911T053324299719Z-fdaa25c7`: DB registration failed before writes because
  a source-snapshot-relative socket path was invalid. The successful rerun uses
  the explicitly configured original project DB socket.
