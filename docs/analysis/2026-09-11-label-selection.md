# Public pIC50 availability and screening selection

## Finding

CYP2D6 has a pronounced activity-dependent label-availability boundary. The
training label population and test acquisition population follow different
selection rules. This is a credible explanation to investigate for poor transfer,
not a measured causal attribution of the leaderboard error or an estimated offset.

## Official evidence

The [assay report](https://openadmet.github.io/octant-cyp-inhib-blog-post/), section
“Tiering strategy”, describes a 50 µM primary screen under TDI conditions, followed
by dual-arm 12-point dose responses on high-inhibition hits. The saved official
announcement additionally specifies log-fold-change and FDR hit criteria. Neither
source establishes a single numerical selection cutoff for all public rows.

Test acquisition used chemisimilars of top hits from CYP1A2, CYP2C9 and CYP3A4;
CYP2D6 used that same compound set. Thus its test set was not selected for 2D6
activity. The announcement explicitly calls the DRC training matrix sparse and
the test matrix dense. Source snapshots:

- `references/official/cyp/assay-dataset/20260907T135802Z/snapshot.md`
- `references/official/cyp/announcement/20260907T135802Z/snapshot.md`

The current assay web page was checked on 2026-09-11 and confirms this design.
This information was already in the saved sources; it should have been inspected
earlier, before advancing a missing-label hypothesis.

## Observable data

Use standardized compound keys, public labels and the existing primary-screen
aggregates. Validate those aggregates directly against the original CSV using
Molecule_Name and SMILES mapped through source aliases. All values agree within
1e-12. Public primary data cover 4,376 compounds per enzyme; none are test compounds.
This is not the complete approximately 11,000-compound experimental library.
“Unlabeled” means no public direct-pIC50 value for that endpoint, not proven assay
failure, proven non-selection, or inactivity.

| Enzyme | Public pIC50 count | Also in primary data | AUC of −log2fc for label presence |
|---|---:|---:|---:|
| 1A2 | 1,412 | 1,412 | 0.8961 |
| 2C9 | 1,285 | 1,285 | 0.5958 |
| 2D6 | 1,493 | 1,493 | 0.9998 |
| 3A4 | 2,335 | 1,805 | 0.4114 |

AUC is a descriptive rank statistic for observed availability, not a trained
activity classifier, a causal effect, or model generalization performance.
The 3A4 result is a warning against extending the simple 2D6 selection story to
all endpoints. Its 530 labels without public primary data and strong responses
without labels show that the public datasets do not expose the full lineage.
We have not resolved those individual non-overlaps or all QC decisions.

### CYP2D6 detail

| Primary TDI log2fc | Compounds | Public pIC50 count |
|---|---:|---:|
| < −2 | 631 | 631 |
| [−2, −1) | 567 | 567 |
| [−1, −0.5) | 355 | 295 |
| ≥ −0.5 | 2,823 | 0 |

Labeled primary median −1.7379; unlabeled primary median +0.3763. All labeled
2D6 rows have FDR < 0.05, but 66.2% of unlabeled rows do too: statistical significance
alone does not identify inhibition because it also includes positive responses.
The ranges overlap (labeled maximum −0.5917, unlabeled minimum −0.8728), so this
is not an exact universal cutoff. Labeled direct-pIC50 median is 4.7260, and 129
labeled compounds are below 4 despite primary TDI selection. Do not substitute a
hard direct-pIC50 label from either presence or primary response.

### Chemical coverage

Median maximum Morgan radius-2/2048-bit Tanimoto from test to endpoint-labeled
training compounds is 0.5172/0.5179/0.4711/0.5377 for 1A2/2C9/2D6/3A4.
2D6 has the lowest coverage by this measure, though pool sizes differ and this
is not a matched-size comparison or proof of a causal error source. Corresponding
unlabeled-to-labeled medians are 0.3279/0.3288/0.3291/0.3418.

For 2D6, labeled/test median MW is 344.42/337.36, logP 3.089/2.859 and TPSA
61.02/70.67. These are covariate differences, not evidence of hidden test potency.
Detailed counts, quantiles and molecular nearest-neighbor keys are saved.

## Modeling implication and next experiment

Existing UMAP OOF evaluates only compounds with pIC50 labels. Clustering their
structures does not restore missing low-primary-response label support. It can
therefore fail to measure precisely the extrapolation required by 2D6 test
acquisition. Importance weighting cannot create target observations in that gap.

The next bounded modeling test should evaluate direct-pIC50 error by **fold-safe
predicted primary response**, then test a conservative continuous low-response
correction with sensitivity analysis. Any assumption about missing direct potency
must be explicit; existing selected-label CV alone cannot validate that region.
A hit-expansion-like structural evaluation among observed labels can complement
UMAP, but still cannot supply unknown inactive outcomes. Recovering compatible
measured weak-activity targets would be stronger evidence than hard pseudo-labels.
This investigation does not yet fit such a correction, reopen external augmentation,
train activity classifiers, or submit predictions.

## Run records and outputs

- Schema probe: `20260911T094230194151Z-2bfa1fd7`.
- Availability, response and chemistry audit: `20260911T094332825500Z-375a66b6`.
  Outputs: availability.parquet, coverage.csv, associations.csv, primary-bins.csv,
  distributions.csv, nearest-labeled.parquet, similarity-summary.csv.
- Independent raw aggregation and FDR audit: `20260911T094442931366Z-2b69ec93`.
  Outputs: raw-audit.csv, primary-availability.png, summary.json.

Runs and their structured manifests are stored in the experiment DB. Every
molecular output uses compound IDs and explicit endpoints/enzyme names. The plot
was visually checked; histograms restrict display to [−5,2], tables retain all rows.
`pixi run check` passes. No infrastructure or model changes were made.
