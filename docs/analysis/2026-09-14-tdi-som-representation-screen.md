# TDI site-of-metabolism and representation screen

## Question

For TDI classification, does a representation trained to identify metabolic
reactivity add information beyond a parent-molecule Morgan fingerprint? The
screen distinguishes atom-level site-of-metabolism (SoM) models, enumerated CYP
products, and a ligand-only machine-learned potential (MLP) embedding. None of
these is treated as a simulation of Compound I, charge/spin transfer, or a full
bond-breaking reaction path.

All downstream comparisons use the frozen UMAP and scaffold folds. Calls use
the leaderboard-inferred rates of 6.533% for CYP2D6 and 33.2% for CYP3A4, with
target prevalences 7.063% and 29.16%. No blind labels are used.

## Sources and leakage controls

- SMARTCyp 2.4.2 was pinned at MDStudio_SMARTCyp commit `92ebd48` and JAR
  SHA256 `d39456ab`. Run `5e7f9e02` produced atom-level generic/3A4 and 2D6
  scores for 6,891 of 6,896 compounds; five charged structures remain missing.
- ATTNSOM was pinned at commit `0faa38b`. Audit `568358Z-9a1e5830` found 98
  challenge connectivities in its nine source SDF files, all on challenge train
  and none on blind test. Full training `221787f2` excluded every one of the
  6,889 challenge connectivities before fitting five folds. It retained 1,679
  molecule-isoform instances over 580 unique external molecules. Its internal
  SoM top-1 accuracy is .747 and top-3 accuracy is .915. Run `c5c2c254`
  extracted 2D6/3A4 atom probabilities and 512-dimensional mean/max-pooled
  FiLM representations for all 6,896 challenge compounds.
- MACE-OFF23 small was pinned at model commit `91a78c5` and SHA256
  `165cce4c`. Run `e80256e2` extracted a 384-dimensional mean/max pooled atomic
  representation from one deterministic RDKit ETKDGv3 conformer for 6,894
  compounds. Two conformer failures remain missing. This is a neutral,
  closed-shell, ligand-only physical representation, not a CYP or TDI model.
- BioTransformer 3.0 was pinned at runnable-JAR commit `6432cf8`, with JAR
  SHA256 `c70cad91`. Probe `be1a4356` generated six annotated one-step CYP
  products for thymol; ten-molecule batch `6d00eba7` generated 103 products and
  mapped every product back to its parent. The full product screen is recorded
  separately below.
- GLORYxR was pinned at commit `ec8ca2b`. Audit `002c2827` found 26 ADAMS rules
  plus three 260-row rule tables, but its only bundled model is explicitly a
  rudimentary tutorial model. The notebook still marks the full FAME3R Zenodo
  model archive as “to be published,” so no scientific model comparison was
  run.
- CYP-MAP was pinned at commit `e55ddd9`; the five released trial checkpoints
  have distinct registered hashes. Source audit `311f1de4` found 155 challenge
  connectivities among 1,119 unique source structures, all on challenge train
  and none on blind test. Frozen extraction `17d2a645` processed all 6,896
  compounds without failure. It records 88 ensemble substrate/atom-SoM/
  bond-SoM/reaction-type summaries and the 512-dimensional molecular input to
  CYP-MAP's nine-CYP substrate head.

## Controlled results

The common LightGBM readout in run `ca43af17` uses identical hyperparameters
for Morgan, SMARTCyp, ATTNSOM, and their concatenations.

| Endpoint | Split | Representation | AUROC | AP | target-prior MCC |
|---|---|---|---:|---:|---:|
| 2D6 | UMAP | Morgan | .5602 | .2715 | .0428 |
| 2D6 | UMAP | Morgan + SMARTCyp | **.5796** | **.3070** | **.0849** |
| 2D6 | UMAP | Morgan + SMARTCyp + ATTNSOM | .5536 | .2605 | .0428 |
| 2D6 | scaffold | Morgan | .5935 | .2905 | .1015 |
| 2D6 | scaffold | Morgan + SMARTCyp | **.6360** | **.3638** | **.1327** |
| 2D6 | scaffold | Morgan + SMARTCyp + ATTNSOM | .5950 | .2930 | .0630 |
| 3A4 | UMAP | Morgan | .6580 | .3149 | .1879 |
| 3A4 | UMAP | Morgan + ATTNSOM | .6956 | .3827 | .2537 |
| 3A4 | UMAP | Morgan + SMARTCyp | .7154 | .3755 | .2786 |
| 3A4 | UMAP | Morgan + SMARTCyp + ATTNSOM | **.7378** | **.4196** | **.3035** |
| 3A4 | scaffold | Morgan | .7307 | .3825 | .3143 |
| 3A4 | scaffold | Morgan + ATTNSOM | .7510 | .4196 | .3376 |
| 3A4 | scaffold | Morgan + SMARTCyp | .7525 | .4329 | .3345 |
| 3A4 | scaffold | Morgan + SMARTCyp + ATTNSOM | **.7758** | **.4617** | **.3701** |

ATTNSOM's scalar max/top-k SoM summaries are weak; its learned representation
is the useful transfer object. For 3A4, SMARTCyp and ATTNSOM are complementary
on both splits. For 2D6, ATTNSOM degrades the SMARTCyp result and should not be
included merely for symmetry.

MACE-OFF readout `7b1b64c2` is weaker and less stable. Morgan plus MACE raises
3A4 UMAP AUROC/AP from .6580/.3149 to .6777/.3704, but scaffold AUROC is nearly
flat (.7307 to .7329) and target-prior MCC falls .3143 to .3035. The same block
is essentially neutral for 2D6. Raw MACE energy and atom-count properties are
not useful. This supports using an MLP embedding as a diversity feature, but
not regarding it as a direct reaction model.

## BioTransformer product screen

Full one-step combined-mode run `abb7661e` generates products for all 6,896
parents. Robust aggregation `9e3ff5cd` retains 69,615 unique parent-product
pairs. Eighteen of 72,418 raw records name an intermediate metabolite rather
than an original challenge molecule as precursor and cannot be assigned
uniquely; they are explicitly omitted. The prespecified features are
product/reaction counts, 2D6- and 3A4-attributed product counts, mass-change
summaries, coarse reaction-name counts, and the mean Morgan vector over unique
predicted metabolites.

Readout `4ba0c8d0` rejects the block. For A4, Morgan plus all BioTransformer
features moves UMAP AUROC/MCC .6580/.1879 to .6649/.2052, but scaffold moves
.7307/.3143 to .7309/.2910. For D6 it moves UMAP .5602/.0428 to
.5746/.0125 and scaffold .5935/.1015 to .6104/.0677. Product features have
some standalone structural signal, but no stable incremental TDI value across
the two splits.

## CYP-MAP follow-up

CYP-MAP is especially relevant to the original question because it has an
explicit molecule-level substrate head for each CYP in addition to atom/bond
SoM and reaction-type heads. The isolated four-feature 2D6/3A4 substrate-score
block is weak, however: direct CYP substrate status is not a proxy for TDI.
The transferable object is again the learned representation.

Readout `7d2bfced` finds that Morgan plus the CYP-MAP embedding raises A4
target-prior MCC from .1879 to .3206 on UMAP and from .3143 to .3685 on
scaffold. The corresponding AUROC/AP values are .7391/.4323 and .7672/.4668.
It does not provide a stable D6 increment.

Complementarity readout `c0dc3f32` adds CYP-MAP to the prior SMARTCyp+ATTNSOM
block. On A4, the embedding-only addition reaches AUROC/AP/MCC
.7501/.4494/.3469 on UMAP and .7814/.4758/.3855 on scaffold. Adding the scalar
substrate/SoM/reaction summaries as well gives .7483/.4409/.3453 and
.7836/.4771/.3932. Conditional grouped bootstrap `7f96bff8` gives a positive
UMAP MCC-delta interval of [.0132,.0702] for the all-feature addition; the
scaffold interval [-.0074,.0542] still crosses zero. This is a strong new A4
model-bank result, but the interval excludes representation selection and
refitting uncertainty.

## Decision

The screen establishes a real mechanism-adjacent representation gain for TDI,
especially CYP3A4. Retain endpoint-specific banks:

1. CYP2D6: Morgan + SMARTCyp is the cleanest tested combination.
2. CYP3A4: add frozen CYP-MAP to Morgan + SMARTCyp + challenge-disjoint ATTNSOM;
   it improves both point estimates and has a positive conditional UMAP-group
   interval, while scaffold uncertainty still crosses zero.
3. MACE-OFF remains a lower-priority diversity block; do not call it a CYP
   reaction model.
4. Reject BioTransformer product features for the current bank.
5. Do not use GLORYxR's tutorial model as if it were the published full model.

These are representation-screen results, not evidence that either combination
automatically replaces the current live Caruana ensemble. The fixed-weight
follow-up `ff991b3e` finds a small possible complement: a 3:1
Caruana/SMARTCyp-ATTNSOM rank blend raises target-prior MCC .37470 to .38087,
AUROC .77571 to .78094, and AP .46969 to .48307. It swaps 11 of 249 blind-test
calls in each direction. The grouped conditional interval for the MCC change is
[-.01363, .02553], so this is a retained diversity probe rather than sufficient
evidence for a candidate or submission.

The newly discovered CYPMol and MetaCYP releases are higher-priority future
checks than additional generic MLPs. CYPMol directly models CYP substrates,
inhibitors and BoMs; MetaCYP fuses CYP ESM-2 and Uni-Mol bond representations.
Both publish code and weights, but were identified after the prespecified
CYP-MAP screen and remain unevaluated here.

## Current-live mechanism blend

At the user's direction, corrected run `20260914T152216279858Z-120542d5`
tests a fixed 25% mechanism-rank addition to the currently scored endpoint
ranks. Failed precursors `20260914T151613162734Z-3732d234` and
`20260914T151636713388Z-c821f451` preserve corrections to snapshotted input
names and the source of UMAP groups. Run
`20260914T151707304794Z-428271d1` is also rejected: it mistakenly used the
blind-test totals of 49 and 249 calls across the entire, larger OOF cohorts,
rather than applying the rates 49/750 and 249/750 within each fold.

With the correct rates, neither addition improves the current model. For
CYP2D6, adding Morgan+SMARTCyp to the deployed ComplementNB rank changes
target-prior MCC `.09988→.08987`, although AUROC `.57926→.58566` and AP
`.29124→.31657` improve. The blind-test vector swaps 13 of 49 calls in each
direction. The grouped MCC-delta interval is `[-.04046,.02311]`, with only
25.7% improvement probability.

For CYP3A4, adding the expanded SMARTCyp+ATTNSOM+CYP-MAP-all bank to the
currently deployed Caruana/PXR-Admetica rank changes target-prior MCC
`.45467→.44240` and AUROC `.80687→.80556`, despite AP rising
`.51391→.52002`. Its interval is `[-.02929,.00498]`, with 8.3% improvement
probability. The mechanism bank's earlier positive result was relative to
Morgan or the pre-PXR Caruana rank; it is not incremental to the current
PXR/Admetica candidate.

Packaging run `20260914T151809559565Z-63dc308a`, candidate
`20260914T151824703191Z-0c3c0267`, and prepared submission
`20260914T151828717456Z-361f4c44` were created from the invalid-rate
interpretation. The file passed schema validation but was never uploaded and
is withdrawn. Retain the currently scored classification vector unchanged.

## Metabolite activity and embedding follow-up

The original BioTransformer screen treated metabolites as an unweighted bag of
fingerprints. Follow-up run `20260915T065421221917Z-7d41918f` instead fits 40
outer-fold-safe Morgan regressors for direct pIC50 and primary log2fc across the
two endpoints and two split strategies. Each fold-specific transform excludes
the entire downstream held-out fold before scoring the parent and its predicted
products. It scores 73,954 unique structures and aggregates 100 features over
69,615 parent-product pairs. The earlier run
`20260915T065250424452Z-2641826c` is superseded because it emitted only held-out
rows rather than a complete fold-routed feature table.

Readout `20260915T065610180920Z-f532943e` finds a modest A4 signal relative to
Morgan. Morgan plus endpoint-specific product-activity summaries raises
target-prior MCC from `.18788→.23022` on UMAP and `.31435→.33295` on scaffold;
AUROC changes `.65799→.66748` and `.73072→.73449`. D6 MCC does not improve on
both splits. More importantly, fixed live-rank audit
`20260915T065852641674Z-0d7e9a75` rejects every A4 addition. The least harmful
10% activity-rank addition changes current MCC `.45467→.44700`, with grouped
interval `[-.01863,.00359]` and 8.8% improvement probability. Failed precursor
`20260915T065828182617Z-ea5315e1` preserves an input-name correction.

Feature run `20260915T070040196804Z-517e5bfc` then embeds all parents and
products with the pinned frozen 2,048-dimensional CheMeleon encoder, and stores
parent, all-product, D6-product and A4-product means plus parent/product distance
summaries. Fold-local PCA readout `20260915T070440285192Z-321d27f7` shows that
the A4 parent embedding itself is useful: Morgan plus parent PCA32 raises MCC
`.18788→.25053` on UMAP and `.31435→.32210` on scaffold. Adding the
endpoint-product-minus-parent embedding reaches `.23803/.32830`; it therefore
does not consistently improve on the parent embedding. All D6 embedding recipes
worsen MCC. Failed run `20260915T070416092824Z-be0da3d1` preserves a corrected
endpoint-to-scope mapping.

Finally, current-live A4 audit `20260915T070659432575Z-4d3886fb` rejects the
embedding ranks. A 10% parent-embedding addition changes MCC `.45467→.44240`,
with interval `[-.02474,-.00148]` and 1.2% improvement probability. The
parent-plus-product-delta form reaches the same point MCC, interval
`[-.02271,-.00228]`, and 0.83% improvement probability. Thus both activity
scoring and embeddings rescue part of a weak Morgan baseline, but the predicted
metabolite contribution is not incremental to the current A4 ensemble. Do not
create a candidate or submission from this follow-up.

Signal diagnostic `20260915T071511594577Z-e24e3cfe` adds outer-safe target
associations, normalized LightGBM gain/split importance and held-out block
permutation importance. The largest individual activity-feature correlations
are small: absolute Spearman rho is about `.09-.10` for both endpoints. A4 TDI
is positively associated with predicted minimum direct pIC50 over products
(`rho=.092/.092` UMAP/scaffold) and negatively associated with maximum predicted
A4 log2fc (`-.093/-.101`). Its predicted parent A4 direct pIC50 has only
`rho=.061`, despite being the largest individual activity gain feature at about
5.8-6.1% of total model gain. This difference confirms that the tree uses
nonlinear thresholds and interactions rather than only marginal correlation.

For the A4 all-activity model, Morgan receives 34-35% of normalized gain. The
remaining gain is distributed across redundant activity blocks; held-out
permutation identifies Morgan as dominant (AUROC drop `.081/.128`) and the
largest metabolite block, all-product absolute scores, as much smaller
(`.018/.012`). The full parent embedding is genuinely informative: it receives
about 71% gain and its permutation lowers A4 AUROC by `.134/.156`. The
A4-product-minus-parent embedding also receives about 60% gain, but its drop is
only `.063/.085`, below both the parent embedding and Morgan in that readout.
Simple product count/distance scalars have weak marginal correlations
(`|rho|≤.053`). D6 shows some AUROC importance but no stable fixed-rate MCC gain,
consistent with an unstable operating-point signal rather than a useful call
vector. Run `20260915T071335176484Z-44b9ddf7` is superseded because its first
gain table averaged columns within blocks instead of summing each block within
fold; correlations and permutations were unaffected.
