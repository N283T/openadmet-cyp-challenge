# Broader ADME pretraining: verified scope and modest next additions

Our current four-source corpus is CYP-focused. It has not yet tested the broader
ADME supervision described for Jeremy's `chemprop_medium`. The most direct next
addition is **five explicitly identified AstraZeneca physicochemical/ADME assays
already present in the local ChEMBL mirror**, followed by the public Biogen
six-endpoint panel. These are proposals, not claims about Jeremy's exact corpus.

## What Jeremy's public artifacts actually establish

The local public repository is
`/home/nagaet/ghq/github.com/jeremycheminf/openadmet_scripts`, revision
`d3a5549e0c98167435cbf8f5115c68a9a184ed7a`.
Its [CYP README](https://github.com/jeremycheminf/openadmet_scripts/tree/main/CYP_Challenge)
states that the sibling pretraining project used public ChEMBL, TDC and Polaris
data plus a public Novartis-surrogate CYP panel. It does not enumerate the
pretraining task names, molecule list or per-source extraction rules. PubChem-only
adaptation is described as a subsequent separate experiment, not the definition
of `chemprop_medium`.

Recorded checkpoint inspection:

- Initial safe inspection `20260912T054633899789Z-1ddb4f66` records the pickle
  global types and strings when the default weights-only loader rejects the
  Lightning AttributeDict. No unsafe loader fallback was used.
- **`20260912T054715454522Z-132eb4d1`** loads using explicitly allowlisted
  installed Chemprop/Torch classes and maps AttributeDict to a plain dictionary.
  It saves numeric architecture metadata and tensor shapes without inference.

| Verified checkpoint property | chemprop_medium | chemprop_chemeleon |
|---|---:|---:|
| Message-passing hidden width | 600 | 2,048 |
| Message-passing depth | 4 | 6 |
| Aggregation | NormAggregation, norm=100 | MeanAggregation |
| Predictor hidden width | 600 | 1,024 |
| Predictor hidden layers | 3 | 3 |
| Predictor task count | 47 | 47 |
| Batch normalization | false | false |
| output_columns | null | null |

The **47 output heads are verified; their names are not recoverable from the
provided metadata**. In particular, “47 heads” must not be converted into an
invented list of 47 datasets. Our depth-5 mean-aggregation scratch encoder is also
not an exact architectural reproduction of his medium checkpoint.

## First addition: five AstraZeneca assays already available locally

Read-only database inventory **`20260912T054818587463Z-14c806d1`** queries
ChEMBL_36 source **27**, whose source record describes AstraZeneca public DMPK
and physicochemical measurements. SQL and all grouped counts/descriptions are
preserved in the run. No new molecule-level training data were exported.

| Assay | Property | Rows / ChEMBL molecules | Important units/conditions |
|---|---|---:|---|
| CHEMBL3301363 | LogD7.4 | 4,200 / 4,196 | Octanol/water, pH7.4, shake flask; already logarithmic |
| CHEMBL3301364 | Solubility | 1,763 / 1,763 | Solid starting material, pH7.4; standardized nM |
| CHEMBL3301365 | Human plasma protein binding | 1,614 / 1,614 | Percent bound, equilibrium dialysis |
| CHEMBL3301370 | Human liver microsomal intrinsic clearance | 1,102 / 1,102 | mL/min/g; 274 lower-censored and 84 upper-censored |
| CHEMBL3301372 | Human hepatocyte intrinsic clearance | 408 / 408 | uL/min/10^6 cells; 104 lower-censored and 15 upper-censored |

This first bundle contains **9,087 observation rows before structure curation and
challenge-overlap exclusion**, not 9,087 unique new compounds. LogD already has
four repeated molecule records; exact assay/molecule replicates need an explicit
aggregation policy. Other assays share compounds with each other and may overlap
the current external corpus.

Practical curation points:

- These assay metadata have confidence score **0**, because single-protein target
  confidence is inapplicable. LogD and solubility have assay type **P**. Reusing
  our CYP single-protein/confidence>=8 policy would wrongly remove them. Use an
  explicit endpoint/assay allowlist and retain source assay conditions.
- Keep logD as reported. A proposed solubility transformation is log10 molar
  solubility from positive nM values, preserving the original concentration and
  assay pH. This is **not** inhibitory pIC50.
- Keep human microsome and human hepatocyte clearance as separate targets;
  retain censor relations and use interval/ranking supervision where appropriate.
  mL/min/g is numerically equivalent to uL/min/mg, while per-million-cell units
  remain a different measurement basis.
- Several source-27 uncensored entries have a null standard relation rather than
  `=`. Treating those as exact needs an explicit deposition-specific policy,
  supported by their source assay description; never silently reuse the IC50
  equality filter.
- Human PPB percent **bound** is not the same field as Biogen log percent
  **unbound**. Preserve the distinction and use different source/measurement heads.

The local mirror has **612 rat hepatocyte** observations in CHEMBL3301371,
separate from the 408 human observations. The TDC page advertises a 1,020-row
hepatocyte dataset; the equality 408+612=1,020 suggests species aggregation and
should be checked at import, not labeled entirely human by assumption.
[TDC source catalog](https://tdcommons.ai/single_pred_tasks/adme/#clearance-astrazeneca).
Do not import TDC wrappers of this AZ deposition as additional independent data.

## Second addition: Biogen's public six-endpoint panel

The authors provide **3,521 commercially sourced compounds** with six measured
endpoints: human/rat microsomal stability, MDR1-MDCK efflux ratio, solubility,
and human/rat plasma protein binding. This supplies diverse auxiliary supervision
in one bounded acquisition. The source CSV contains compound/vendor identifiers
and logged properties. Use the author repository as the primary snapshot and
preserve its commit/hash, header definitions and original missing values.
[Author repository and raw-data description](https://github.com/molecularinformatics/Computational-ADME).

Polaris distributes a curated 3,516-molecule form of this same panel; **it is an
alternative representation, not a second independent dataset**. Its catalog
records log clearance in mL/min/kg, log solubility at pH6.8 in ug/mL, log percent
unbound and log MDR1-MDCK B-A/A-B efflux ratio. These should remain source-specific
heads; MDCK efflux ratio must not be labeled passive permeability.
[Polaris Biogen dataset](https://polarishub.io/datasets/biogen/adme-fang-v1).

The first experiment can add these six heads with a fixed total ADME-family
sampling budget, keeping CYP source exposure comparable. All challenge
connectivity blocks should again be excluded before external splitting. Predictive
benefit is a hypothesis: adding unrelated endpoints is not itself evidence that
the CYP representation improved.

## Small optional follow-up, not an initial bulk collection

If an additional CYP mechanism signal is useful, TDC lists Carbon-Mangels
substrate/non-substrate labels for CYP2C9/2D6/3A4 (666/664/667 entries). These are
**substrate classification**, not inhibition potency, and derive from a literature
aggregation; shared molecules and publication provenance need checking. They may
add target-specific supervision but little new chemistry after overlap removal.
[TDC substrate source definitions](https://tdcommons.ai/single_pred_tasks/adme/#cyp2c9-substrate-carbon-mangels-et-al),
[original classification study](https://onlinelibrary.wiley.com/doi/10.1002/minf.201100069).

A broad AqSolDB/Caco-2 collection is lower priority than the two identified
bundles above: we already have a specific solubility addition and a defined
transport readout, and the current comparison needs to isolate whether adding
non-CYP measured properties helps. New tasks should be added in a separate
versioned corpus and controlled training arm.

## What exists in the previous PXR workspace

The earlier PXR project retains the full ChEMBL36 PostgreSQL mirror and archive.
Its recorded broad pairwise preparation report was inspected in run `1ddb4f66`:
250,000 sampled raw bioactivity rows yielded 75,608 activities across 2,900 assays
and 841,412 ordered pairs. Its scope was B/F IC50/Ki/Kd/EC50/AC50 with confidence
>=8, **not** a curated ADME property collection. It cannot substitute directly for
the source-27 physicochemical/clearance extraction above.

A file named `data/admet_ai_predictions.parquet` also exists. It is explicitly
model predictions, not a new set of measured ADME labels. The inventory did not
find an already prepared measured Biogen/AZ broad-ADME training table in the
reviewed project data paths. This is a bounded file/source inspection, not a
claim that no such file exists elsewhere on the machine.

No new training, GPU job or large data download was launched for this inventory.

## Implemented AZ addition after inventory approval

The subsequent requested implementation is complete:

- Raw measured-data extraction: **`20260912T055219774452Z-fa008fc6`**,
  9,087 observations across 5,599 ChEMBL molecules.
- **Curated/merged corpus: `20260912T055348682784Z-be731fe3`**.
  Use `external-compounds.parquet` and `observations.parquet` in its outputs.
- Independent transformation/exclusion audit:
  **`20260912T055456975084Z-4a4a17ab`**, passed.

After standardization and excluding 331 challenge-overlap observation rows,
AZ supplies **5,423 compounds, 8,756 observations, 8,295 exact continuous labels
and 461 censored labels**. Of those compounds, **4,812 are new** relative to the
frozen three-source CYP corpus `20260912T051513413487Z-ac9df6b5`.
The merged corpus contains **40,132 compounds, 165,128 observations, 76,119
continuous labels and the unchanged 125,063 binary labels**.

| Preserved ADME head | Exact continuous | Censored |
|---|---:|---:|
| ADME_LOGD74 | 4,053 | 0 |
| ADME_SOLUBILITY74 | 1,721 | 0 |
| ADME_PPB_HUMAN | 1,539 | 0 |
| ADME_HLM_CLINT | 720 | 347 |
| ADME_HEPATOCYTE_CLINT | 262 | 114 |

All five heads support continuous/ranking learning. No arbitrary classification
thresholds were added. Source `chembl36_az_adme` remains a separate sampling group;
the `isoform` field contains explicit `ADME_*` task categories for compatibility
with the current encoder trainer, not an assertion that these are CYP isoforms.

Transformations are identity for logD, `log10(nM)-9` for molar log-solubility,
`log10(100-percent_bound)` for log percent unbound, and log10 of the separately
recorded native clearance units. PPB transformation is decreasing and reverses
inequalities; the others preserve them. Clearance species are human only.
Null source relations are explicitly treated as uncensored for this five-assay
allowlist; raw relations, assay measurement ranges and the assumption flag remain
in `az-curated-records.parquet`. No additional range-based censoring is invented.

`az-challenge-overlap-observations.parquet` retains excluded overlaps separately.
`az-structure-provenance.parquet`, `head-metadata.json`, the raw acquisition's
`source.json`, original activity/assay IDs and full curated source records preserve
identity, species, units, source descriptions and transform provenance.

This adds measured non-CYP supervision to an immutable corpus. No encoder training
or GPU work was performed by the data-preparation task, and prior comparison
corpora were not changed.

## Implemented Biogen addition after AZ

The primary-author CSV, README and license were acquired in
**`20260912T055605172480Z-f639f5a2`**, resolved to commit
`b00df003de117ce9e5b381afd886095c5f2af2d5` of
`molecularinformatics/Computational-ADME`. CSV SHA-256:
`2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b`.
The raw file is 388,279 bytes and contains 3,521 compounds and 11,318 nonmissing
property measurements. No Polaris copy was added.

Schema/censor inspection `20260912T055622514051Z-3be11444` finds **no explicit
inequality tokens** in the six property columns. Values are already logged and
are retained numerically unchanged; repeated plateau values are not converted
into invented censor bounds. The primary raw snapshot has only 194 human and
168 rat PPB measurements before overlap exclusion, fewer than other endpoints.

**Merged AZ+CYP+Biogen corpus: `20260912T055809060996Z-06088efb`.**
Use its `external-compounds.parquet` and `observations.parquet`.
Biogen adds 3,483 usable compounds, of which **3,364 are new**, and **11,208
continuous observations** after excluding 110 challenge-overlap observation rows.
All six heads support continuous/ranking learning; no binary thresholds are
fabricated.

| Biogen head | Accepted exact observations | Original property units |
|---|---:|---|
| Human microsomal clearance | 3,062 | Logged mL/min/kg |
| Rat microsomal clearance | 3,029 | Logged mL/min/kg |
| MDR1-MDCK efflux ratio | 2,612 | Logged B-A/A-B ratio |
| Solubility pH6.8 | 2,156 | Logged ug/mL |
| Human plasma protein binding | 189 | Logged percent unbound |
| Rat plasma protein binding | 160 | Logged percent unbound |

The merged corpus has **43,496 compounds, 176,336 observations, 87,327 exact
continuous labels and unchanged 125,063 binary labels**. Biogen assay names,
organism where specified by the endpoint, property units, original column/token,
identity transform, source commit and CSV hash are retained. MDR1-MDCK has an
explicit assay-system field and unspecified organism rather than an invented
species assignment. It remains distinct from passive permeability and from
microsomal clearance.

Raw/curated observations and excluded overlaps are saved separately. Independent
verification **`20260912T055905891725Z-f1e5161d`** confirms every one of the
11,208 accepted values exactly equals its original numeric logged token, verifies
the input CSV hash and the six unique keyed heads, and confirms no challenge
connectivity or fabricated binary labels. This data-preparation task used no GPU
and did not modify earlier frozen corpora or trainer behavior.

## Frozen upstream CheMeleon control

Export run **`20260912T061231617816Z-299439e8`** extracts all 6,896 challenge
compounds from Jeremy's upstream `checkpoints/chemprop_chemeleon.pt`. The safe
`weights_only=True` loader and strict message-passing state load verify 2,048
hidden dimensions, depth 6, MeanAggregation, no batch normalization and **47
pretraining output heads**. The upstream README and script03 are registered
inputs: they identify this checkpoint as the sibling ADME-pretraining artifact,
separate from the challenge's 17-head multitask fine-tune. Original training
molecule manifests and named 47-task definitions remain unavailable, as with the
upstream medium checkpoint. Extraction uses no challenge labels and adds no
fine-tuning.

LGBM run **`20260912T061253921227Z-4930ca7c`**, matched audit
**`20260912T061613446259Z-edd55762`**, and recorded comparison
**`20260912T061638025159Z-a9d3b413`** use the same retained parent CV and all
6,384 OOF label rows. These are plain frozen-embedding LGBM readouts, without
auxiliary log2fc features or PCA, using identical LGBM settings.

| Upstream representation | Overall ST-RAE | Neighbor ST-RAE |
|---|---:|---:|
| CheMeleon 2,048 | 0.674158 | 0.786083 |
| Medium 600 | 0.651179 | 0.739368 |

Medium performs better in this controlled readout. This comparison does not
establish the best TabPFN/TabICL readout or the effect of vanilla CheMeleon
initialization alone: the two upstream checkpoints also differ in architecture
and unpublished pretraining details. Frozen embeddings, all LGBM models, OOF
predictions and matched audit tables are saved in the runs above.

### CheMeleon with nested log2fc and the prespecified four-way mean

Run **`20260912T062323161344Z-db6745c3`** adds the four nested log2fc values
from `bac5fb3a` to the 2,048 frozen upstream CheMeleon features, without PCA,
using the same fixed LGBM settings and retained parent CV. Audit
**`20260912T062641940553Z-33400514`** gives overall ST-RAE **0.645076** and
neighbor ST-RAE **0.739819** (plain CheMeleon was 0.674158 / 0.786083).

The prespecified equal four-way mean combines original Caruana `96fce32b`,
Jeremy medium frozen TabPFN `2f0fbfbd`, Jeremy medium fine-tuned log2fc LGBM
`730c8bea`, and this new CheMeleon model. Blend
**`20260912T062715715731Z-1fd5fd13`**, audit
**`20260912T062716598528Z-b318390d`**, gives **0.605600 overall / 0.696930
neighbors**. No OOF weights were fitted. Overall performance does not improve
on the preceding three-way mean (0.605039).

Matched diagnostic **`20260912T062742470988Z-460d1539`** verifies all 6,384
molecule/endpoint/fold keys and label/confidence/neighbor fields, then records
prediction and residual Pearson correlations:

| CYP | CheM vs medium FT prediction | CheM vs medium FT residual | CheM vs medium frozen PFN residual |
|---|---:|---:|---:|
| 1A2 | 0.8946 | 0.9546 | 0.9215 |
| 2C9 | 0.9439 | 0.9442 | 0.8826 |
| 2D6 | 0.8460 | 0.9783 | 0.9462 |
| 3A4 | 0.9701 | 0.9440 | 0.8490 |

The strong residual correlations are consistent with the limited benefit of
adding CheMeleon to this particular three-model mean. They do not establish that
all readouts of this representation are redundant. This control used CPU only;
no further training was started after the results were available.
