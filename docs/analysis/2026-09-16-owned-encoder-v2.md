# Project-owned encoder v2: Jeremy-medium replacement and complement screen

The second project-owned encoder program closes almost all of the frozen
`chemprop_medium.pt` gap. Under the matched raw-embedding LightGBM protocol,
target-specific owned choices beat Jeremy on CYP1A2, CYP2D6 and CYP3A4. CYP2C9
remains weaker as a standalone replacement, but is a useful independent peer in
the existing target-specific flat prediction bank.

This is a development-CV result. Target-wise selection was made after viewing
the same outer OOF scores, so the selected four-target macro is descriptive and
not an unbiased estimate of future performance.

## Data and leakage policy

Source acquisition runs `20260916T081236163734Z-3ae85fcd` and
`20260916T081236163734Z-aa990e05` added the complete public TDC ADME registry
and the Polaris/Novartis CYP panel. These were combined with the existing local
ChEMBL 36, PubChem, Tox21, Biogen, EPA NVS and OpenADMET source families.

Curation run `20260916T082302554619Z-656c63e6` standardizes the sources into
**78,449 compounds** and removes **4,421 observations** whose standardized
connectivity overlaps the challenge universe. It retains four task-layout
views:

| View | Observations | Meaning |
|---|---:|---|
| assay | 388,673 | Source assay heads retained |
| compact-all | 388,673 | ChEMBL consolidated by isoform and measurement type |
| compact-CYP | 283,994 | CYP-related compact subset |
| endpoint-only | 356,143 | Assay-fragmented ChEMBL excluded |

Three known TDC duplicates of separately acquired AstraZeneca datasets are not
counted as independent training sources. The public source families are now
represented, but Jeremy's exact molecule list, 47 tasks and curation remain
unknown and are not claimed to have been reproduced.

## Matched evaluation

All primary comparisons use the same five parent folds, challenge labels and
official soft-threshold RAE implementation. Encoders are frozen, challenge
molecules are represented by raw graph embeddings, and the same LightGBM
readout is fit independently for each target. Prediction rows carry compound,
endpoint and fold identifiers. PCA and fine-tuning are follow-up ablations, not
mixed into the primary comparison.

The data-layout screen (`20260916T084005484014Z-200a0a24`) reduces the owned
macro ST-RAE from the earlier `.792471` to `.686153` with endpoint-only
supervision. The architecture screen (`20260916T091316402940Z-378864af`)
selects width 600, depth 3, normalized aggregation and a three-hidden-layer
pretraining head at `.681719`. Mean aggregation is consistently poor
(`.725837` to `.733644`).

Raw feature concatenation is the largest remaining gain. Audit
`20260916T091919159581Z-99c5e238` reaches `.656595` with `broad5`, which joins
three endpoint-only architecture views with compact-all and assay views at the
feature level. This is not a prediction ensemble: one LightGBM sees one joined
raw feature matrix keyed by compound.

Readout audit `20260916T095245039121Z-5949843f` compares default LightGBM, five
capacity/regularization variants and fold-local PCA256/PCA512 on `top_arch` and
`broad5`. Raw `broad5_leaves31` reaches macro `.653291`, only `.002112` behind
Jeremy `.651179`. PCA is decisively harmful: the four PCA macros range from
`.707214` to `.723308`.

## Target-specific frozen results

Lower ST-RAE is better. These are the selected owned raw readouts, compared
with the same Jeremy raw600 LightGBM control.

| Target | Selected owned view/readout | Owned | Jeremy | Delta |
|---|---|---:|---:|---:|
| CYP1A2 | broad5, leaves31 | **.720874** | .728514 | -.007640 |
| CYP2C9 | broad5 + 2C9 specialist, default | .552196 | **.535216** | +.016980 |
| CYP2D6 | top-architecture concat, slow31 | **.850983** | .858244 | -.007261 |
| CYP3A4 | broad5, leaves31 | **.474445** | .482742 | -.008297 |

The selected owned-only descriptive macro is `.649624`, versus `.651179` for
one common Jeremy raw600 readout. The meaningful conclusion is target-wise:
the owned representation is already a plausible Jeremy replacement for three
targets, while CYP2C9 is not.

## CYP2C9 specialist and flat-bank placement

Subset run `20260916T094225450353Z-44e1d7d6` builds 59,255 direct CYP2C9
observations across 33,281 compounds and seven consolidated heads from ChEMBL,
PubChem, Tox21, TDC and Novartis. Its width600/depth3 specialist encoder
(`20260916T094226751324Z-20c2d3d6`) is weak alone at `.655998`. Appending its
raw embedding to `broad5` improves the default readout `.563412→.552196` in
matched audit `20260916T095531264237Z-00304122`, but does not beat Jeremy.

The specialist is nevertheless complementary. Flat-bank run
`20260916T095939912871Z-fd5036b3` adds the complete owned readout as one atomic
peer to the CYP2C9-only bank and reruns equal-all, median, held-fold equal
subset, Ridge, positive linear, MSE/interval simplex and Caruana. Equal-subset
improves current control `.463103→.456116`; all five held-out-fold searches
choose the identical four-member set:

- Jeremy scores + log2FC TabPFN
- Jeremy fine-tuned + log2FC LightGBM
- CheMeleon PCA256 + log2FC TabICL
- owned broad5 + CYP2C9-specialist raw LightGBM

Caruana `.456237` and interval simplex `.456556` corroborate the gain and give
the owned peer positive weight in every fold. This follows the flat-bank policy:
no nested prediction member and no 10/25/50% defensive addition.

## Fine-tuning and negative results

Fold-safe interval fine-tuning of the best endpoint-only depth-3 checkpoint is
recorded in `20260916T095627181174Z-7b499414`. Its direct predictions score
`.751903/.560206/.870445/.480192`; the fold-specific LightGBM embeddings score
`.727241/.561548/.866827/.491350` for 1A2/2C9/2D6/3A4. Neither displaces the
best frozen target-specific readout. The run is retained because it shows that
the frozen external representation, not challenge-label adaptation, carries
the present gain.

Other informative negatives are:

- compact-CYP supervision is worse than broader task layouts (`.718911` macro);
- width 384, 768 and 1024 do not beat width600/depth3;
- mean aggregation is substantially worse than normalized aggregation;
- PCA256/PCA512 lose heavily on both leading concatenated views;
- the CYP2C9 specialist is useful only as a complementary view, not alone;
- specialist slow31 (`20260916T095557213469Z-d7f5c308`) is worse than its
  default readout;
- initial specialist audit `20260916T095436053347Z-d8555cf6` failed because
  all-target controls were not filtered to CYP2C9. The retry fixes the generic
  audit path and is the authoritative result.

## Decision

Retain `broad5_leaves31` for CYP1A2 and CYP3A4 and `top_arch_slow31` for CYP2D6
as project-owned frozen replacement candidates. Retain the specialist-augmented
readout as a full CYP2C9 bank peer, not as a standalone Jeremy replacement.
Do not use PCA for these owned embeddings, and do not promote the fine-tuned
variant. No inference, candidate package or leaderboard submission was created.
