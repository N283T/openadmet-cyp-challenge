# Research log

The canonical research narrative lives here, in English, tracked in this private
Git repository. This carries forward the dated research-log format previously
used in PXR issues #100 and #208, without requiring an Issue for each experiment.

## Current direction

- **Checkpoint predictor views produce a stronger D6 flat bank:** six
  distillation/adaptation checkpoints now expose raw600, final-hidden,
  external-task score and raw-plus-score views. Across 32 atomic readouts, the
  59-task anchored raw-plus-score LightGBM is the best D6 single at `.835963`,
  while target winners differ elsewhere. Adding six complete checkpoint peers
  to the established D6 bank improves outer-cross-fit equal subset from the
  current control `.828397` to `.815840`. Keep the other three established
  banks unchanged. Full-fit deployment produced three valid D6-only candidates
  on the favorable scored placement baseline. Prefer stable-five candidate
  `40e0984c` (OOF `.811998`); equal-ten `f1786ea3` is the conservative control
  and checkpoint-only `ed3d1154` is exploratory. None is uploaded. See
  [RL0112](2026-09.md#rl0112--predictor-views-and-flat-checkpoint-ensembles-improve-d6)
  and [RL0113](2026-09.md#rl0113--three-full-fit-d6-checkpoint-submission-candidates).

- **Representation distillation produces the first universal student that
  beats frozen Jeremy macro:** a random teacher-only student reaches `.647986`
  and beats Jeremy on 1A2, 2C9 and 3A4. A Jeremy-initialized, externally
  supervised anchored student improves further to macro `.646348`, driven by
  2C9 `.519131` and 3A4 `.474981`; Jeremy remains best on D6. The owned-init
  hybrid is negative. Keep the two successful students as separate atomic,
  target-specific candidates. A Jeremy-initialized PubChem-only follow-up moves
  D6 only `.858244→.858168`, while teacher anchoring worsens it; retain frozen
  Jeremy for D6. See
  [RL0110](2026-09.md#rl0110--jeremy-representation-distillation-surpasses-the-frozen-teacher)
  and [RL0111](2026-09.md#rl0111--jeremy-initialized-pubchem-only-adaptation-is-marginal-for-d6).

- **Challenge train plus auxiliary adaptation gives a new owned 3A4
  candidate:** fold-safe 12-head adaptation uses direct pIC50, primary log2FC
  and direct-inhibition Emax together. Raw LightGBM scores
  `.737444/.554001/.877706/.480498`; 3A4 beats Jeremy `.482742`, and 2C9
  improves slightly, while 1A2 and D6 worsen. Keep it target-specific and make
  no bank change yet. See
  [RL0109](2026-09.md#rl0109--fold-safe-challenge-train-and-auxiliary-adaptation).

- **Joint sparse-matrix training nearly closes the owned representation gap:**
  a from-scratch 59-task continuous model nearly matches Jeremy on 1A2, while
  a 28-task CYP model nearly matches it on D6. Combining those with the earlier
  target-specific owned winners gives `.729526/.554558/.859557/.488025`,
  descriptive macro `.657917`, versus Jeremy `.651179`. Source layouts chosen
  to resemble Jeremy are worse, and including binary tasks also hurts: use the
  checkpoint analysis as a mechanics clue, not a reproduction target. Keep
  target-specific candidates separate and make no bank change yet. See
  [RL0107](2026-09.md#rl0107--joint-sparse-matrix-training-nearly-closes-the-owned-gap)
  and the [report](../analysis/2026-09-16-jeremy-medium-gap-reanalysis.md#joint-matrix-experiments).

- **Small refinements do not finish the single-checkpoint replacement:**
  matched 2x/4x target-isoform loss weights worsen every intended target.
  Honest nested LightGBM selection scores owned versus Jeremy at
  `.734005/.551885/.860291/.490499` versus
  `.726729/.531426/.856142/.478637`; Jeremy also benefits from tuning. Local
  25/50/75% encoder-weight interpolation produces no new target winner. The
  next attempt should change the pretraining signal or formulation rather than
  continue loss/readout micro-tuning. See
  [RL0108](2026-09.md#rl0108--small-joint-encoder-refinements-do-not-surpass-jeremy).

- **The earlier one-task-sampled encoders remain useful target-specific
  candidates:** the leak-excluded v2 corpus contains 78,449 compounds. Its
  task/source optimization produced target-specific values
  `.736899/.554558/.870970/.488025` still trail Jeremy
  `.728514/.535216/.858244/.482742`. The earlier 2C9 specialist-augmented
  readout remains a complete flat-bank peer, while the 2C9 and 3A4 models remain
  the best project-owned choices after the joint-matrix screen. See
  [RL0101](2026-09.md#rl0101--project-owned-encoder-v2-closes-the-jeremy-medium-gap)
  and [RL0105](2026-09.md#rl0105--full-oof-task-design-audit-and-follow-ups).

- **Regression now has a flat atomic prediction-bank interface:** equal-all,
  median, held-fold equal subset, Ridge, positive linear, MSE/interval simplex,
  and Caruana share one outer-cross-fit protocol. The current aggregate is a
  control, never a member. Selection is endpoint-specific and macro is only a
  descriptive summary. D6 favors swapping PCA for raw, with interval simplex
  `.827988` and Caruana `.827973`; this is not transferred to the other three
  targets. Four separate bank manifests now cover 1A2, 2C9, 2D6 and 3A4; the
  latter includes separately refit class- and Morgan-quantile LightGBMs rather
  than their nested average. Bank additions and family swaps can now be rerun
  target by target through the unchanged auditable suite. See
  [RL0097](2026-09.md#rl0097--regression-chemeleon-raw-readout-and-full-member-placement).

- **Regression raw CheMeleon is a D6-specific full-member candidate:** matched
  raw2048 TabICL improves PCA256 macro ST-RAE `.650170→.637270`, and their
  standalone equal mean reaches `.627101`. Uniformly adding raw to the current
  endpoint ensembles worsens macro `.594440→.595178`, but exhaustive
  held-fold-excluded equal-member selection chooses the same raw-containing
  four-member D6 subset in all five folds and improves `.828397→.825571`.
  Retain it for exploratory D6 follow-up; do not deploy globally. Routine
  ensemble development uses complete members and subset reselection rather than
  10/25/50% additions. See
  [RL0097](2026-09.md#rl0097--regression-chemeleon-raw-readout-and-full-member-placement)
  and the [report](../analysis/2026-09-16-regression-chemeleon-raw.md).

- **Flatten the legacy CYP3A4 nested member:** the scored vector packages a
  weighted two-LightGBM Caruana block as one of five top-level members. Treating
  those LightGBMs as separate peers in a six-member equal mean changes ST-RAE
  only `.398316→.400327` and improves three of five folds; deleting them costs
  more at `.409099`. The two atomic LightGBMs are now present in the 3A4-only
  bank; equal-subset consistently retains Morgan but drops class and scores
  `.400035`.
  The scored vector is unchanged until explicitly repackaged. See the
  [RL0097 follow-up](2026-09.md#rl0097--regression-chemeleon-raw-readout-and-full-member-placement).

- **A4 raw600 helps as a complement, not a wholesale replacement:** matched
  Jeremy600 TabICL/TabPFN raw readouts improve standalone AUROC/AP over PCA256.
  Replacing both PCA members loses ensemble MCC, but retaining both raw and PCA
  yields a seven-member fixed equal mean selected identically in all five outer
  folds. It improves AP/MCC from `.51391/.45467` to `.52275/.46540`, with a
  slight AUROC decrease `.80687→.80578`; its MCC interval still crosses zero.
  Keep it as the first A4 development candidate and a six-member mean as the
  compact fallback. See
  [RL0095](2026-09.md#rl0095--a4-matched-raw-readouts-and-fixed-mean-simplification)
  and the [report](../analysis/2026-09-15-a4-raw-simplification.md).

- **The first D6 candidate is the 11-model raw600/PCA256 block average:** its
  controlled cross-fit AUROC/AP/target-prior MCC is
  `.58076/.29686/.11011`; keep the PCA256-only ten-model variant as the cheaper
  fallback. See
  [RL0094](2026-09.md#rl0094--d6-rawpca-block-prioritization-and-remaining-a4-raw-controls).

- **Classification now has an extensible OOF prediction-bank interface:** a
  prespecified nested comparison evaluates ridge, three convex-simplex losses,
  and a permutation-invariant expert-set network with the current fixed ratio
  retained as an eligible expert/control. No learned D6 stacker beats the
  control MCC `.09988`; the 24-expert equal mean ties MCC and improves
  AUROC/AP. A4 ridge improves all point metrics to `.80924/.53012/.45774`, but
  its grouped MCC interval crosses zero. Anchor-preserving follow-up does not
  rescue D6. Keep the D6 fixed ratio for submission, retain A4 ridge only as a
  development candidate, and add future models as keyed OOF columns rather
  than inventing new manual blend locations. See
  [RL0093](2026-09.md#rl0093--extensible-classification-stacking-and-safe-anchoring).

- **The legacy PCA-only A4 pool compresses to five members:** matched raw
  versus PCA readouts show that raw CheMeleon parent embeddings are strongest
  for A4, while no D6 embedding view transfers across both splits. A
  hierarchy-correct peer audit places SMARTCyp, ATTNSOM, CYP-MAP,
  BioTransformer, raw embeddings and MACE in the same fold-excluded Caruana
  pool; expanded pools worsen or tie both endpoints, and explicit swaps are not
  bootstrap-stable. A4 can drop Morgan balanced and external-class PCA256 from
  the development pool without a material cross-fit loss, but a three-member
  compression loses significantly. RL0095 supersedes this as the active A4
  development direction after adding matched raw members. A corrected genuinely flat D6 follow-up
  performs substantially worse than the hierarchy and does not provide a viable
  successor. Exhaustive held-fold-excluded selection over all 2,047 subsets of
  the 11-member D6 core finds ten total components as the smallest design that
  strictly preserves all three point metrics: AUROC `.57933→.58059`, AP
  `.29059→.29623`, and target-prior MCC `.09988→.11011`. Contrary to the limited
  three-member ablation, its stable core retains Morgan, Mordred and FP+Mordred
  and removes four different legacy members. See
  [RL0092](2026-09.md#rl0092--classification-simplification-raw-embeddings-and-peer-replacement)
  and the [report](../analysis/2026-09-15-classification-simplification-and-replacement.md).

- **MoLFormer pooling plus Random Forest is a stable D6 complement:** complete
  frozen CLS+mean extraction covers all 6,896 compounds. A 25%
  scores+CLS+mean+log2fc4 Random Forest addition improves ST-RAE
  `.822681→.818407` and Spearman `.476425→.482580`, improves all five folds,
  and has a conditional group-bootstrap interval `[-.00859,-.00006]` / 97.6%
  improvement probability. A fixed average with the earlier hybrid CatBoost
  reaches `.816803/.482720`; learned fold-excluded stacks fail every fold.
  The stable standalone Random Forest was full-fit and submitted as
  `20260914T150145548210Z-1ab27b1d`, but the blind result is negative: D6
  ST-RAE worsens `.565520→.572470`, Spearman `.431248→.415591`, and rank
  `21→22`; MA-ST-RAE worsens `.475974→.477712`. Reject this deployment and do
  not infer that the higher-gain post-screen fixed average will transfer.
  KERMT remains nearly flat; defer CYPMol until a
  direct Uni-Mol comparison is justified. Fixed transfer to 1A2 and 2C9 is
  effectively flat (`-.00008` ST-RAE each), while every 3A4 addition worsens;
  treat this as D6-specific rather than a global representation update. See
  [RL0089](2026-09.md#rl0089--d6-molformer-random-forest-leaderboard-submission),
  [RL0088](2026-09.md#rl0088--cross-endpoint-transfer-of-the-d6-selected-molformer-view),
  [RL0087](2026-09.md#rl0087--molformer-pooling-tree-readouts-and-candidate-combinations),
  [RL0086](2026-09.md#rl0086--generic-kermt-and-molformer-cyp2d6-representations),
  and [the report](../analysis/2026-09-14-regression-d6-model-family-screen.md#molformer-pooling-and-candidate-combinations).

- **Ordinal regression is the most stable new D6 classical complement:** a
  second fixed parent-CV screen compares joint interval/pairwise loss, weighted
  pairwise logistic, ordinal logistic, additive splines, Huber, and exact
  Tanimoto kernel ridge. Interval+pairwise25 gives the best pooled ST-RAE
  `.819885`, but improves only three folds. Ordinal-logit5 at 10% reaches
  `.820686` with four improving folds and a conditional group-bootstrap interval
  `[-.00470,+.00070]` / 92.1% improvement probability. Huber corroborates the
  scores51 signal; Tanimoto KRR is diverse but worsens four folds. Retain
  ordinal10 beside XGBoost-pairwise10 as the stable next-candidate shortlist;
  no inference or submission. See
  [RL0085](2026-09.md#rl0085--cyp2d6-ordinal-interval-pairwise-kernel-and-robust-readouts)
  and [the report](../analysis/2026-09-14-regression-d6-model-family-screen.md#ranking-kernel-and-classical-follow-up).

- **Direct pairwise learning adds a modest new D6 ordering signal:** fixed
  parent-CV comparison of five rankers across three views finds that a 10%
  scores51 XGBoost-pairwise addition improves ST-RAE `.822681→.820434` and
  Spearman `.476425→.477486`. Its prediction correlation with latest is only
  `.6738`, although residual correlation remains `.9567`. Scores-only pairwise
  logistic gives more rank gain at 10% (`.479816`) but less loss gain; 25%
  overweights it. The conditional group-bootstrap interval for XGBoost10
  `[-.00534,+.00069]` crosses zero, with 93.1% improvement probability and four
  of five folds improving. Keep compact pairwise rankers as low-weight model-bank
  entries and follow with a prespecified pair-margin/weight comparison plus
  fold-excluded blending; reject the tested RankNet and do not infer or submit.
  See [RL0084](2026-09.md#rl0084--interval-confident-cyp2d6-pairwise-rankers) and
  [the analysis](../analysis/2026-09-14-regression-d6-model-family-screen.md#pairwise-ranking-follow-up).

- **Broad D6 readout screen finds a modest CatBoost complement:** fixed parent-CV
  comparison of 12 classical/modern model families across three frozen feature
  views produces 36 complete OOF vectors. A 50% moment-matched hybrid CatBoost
  addition improves latest D6 ST-RAE `.822681→.817687` and Spearman
  `.476425→.479956`; scores-only Random Forest and PLS give smaller complementary
  gains. The conditional group-bootstrap CatBoost interval
  `[-.01150,+.00097]` crosses zero and residual correlation remains `.9903`, so
  this is a follow-up candidate rather than a new production model. The tested
  compact spline KAN is weak and should only be revisited with a materially
  different formulation. Confirm a small CatBoost family and fold-excluded
  blending next; no inference or submission. See
  [RL0083](2026-09.md#rl0083--broad-cyp2d6-model-family-diversity-screen) and
  [the report](../analysis/2026-09-14-regression-d6-model-family-screen.md).

- **CYP-trained representations provide a real A4 TDI gain:** SMARTCyp atom
  summaries, a challenge-disjoint ATTNSOM embedding, and frozen CYP-MAP are
  complementary. Adding CYP-MAP to Morgan plus the first two reaches AUROC
  .748/.784 and target-prior MCC .345/.393 on UMAP/scaffold, versus Morgan
  .658/.731 and .188/.314. The conditional UMAP-group MCC interval is positive;
  scaffold uncertainty still crosses zero. D6 instead prefers
  Morgan plus SMARTCyp; adding ATTNSOM hurts. A fixed 3:1 live-Caruana/SoM rank
  blend gives a small +.00617 descriptive MCC change, but its grouped interval
  crosses zero. MACE-OFF is mixed, BioTransformer products fail the dual-split
  increment test, and GLORYxR lacks its published full weights. CYPMol is the
  next frozen representation to audit, followed by MetaCYP; aweSOM weights need
  an access request. Retain endpoint-specific representation banks; no candidate
  or submission.
  See [RL0081](2026-09.md#rl0081--site-of-metabolism-and-mlp-representation-screen)
  and [the report](../analysis/2026-09-14-tdi-som-representation-screen.md).

- **Current-live mechanism additions are rejected after a call-rate correction:**
  the first audit incorrectly used the blind-test totals 49/249 across the full
  OOF cohorts. Corrected 49/750 and 249/750 fold rates show that a 25%
  Morgan+SMARTCyp addition lowers D6 target-prior MCC `.09988→.08987`, and the
  expanded CYP-MAP mechanism bank lowers A4 `.45467→.44240`. Earlier A4 gains
  remain valid versus Morgan/pre-PXR Caruana, not the current PXR/Admetica
  vector. Candidate `0c3c0267` / submission `361f4c44` was never uploaded and
  is withdrawn; retain the scored classification vector. See
  [RL0090](2026-09.md#rl0090--d6-smartcyp-live-blend-candidate) and
  [the report](../analysis/2026-09-14-tdi-som-representation-screen.md#current-live-mechanism-blend).

- **Scored BioTransformer metabolites do not improve the live model:**
  fold-safe pIC50/log2fc scores for predicted products and frozen CheMeleon
  parent/product embeddings both improve the weak A4 Morgan baseline on the two
  splits. The gain is mostly attributable to the parent representation, however,
  and fixed 10% additions lower current-live MCC `.45467→.44700` for activity
  features and `.45467→.44240` for embeddings. D6 has no stable gain. Retain
  these as model-bank diagnostics only; create no candidate or submission. See
  [RL0091](2026-09.md#rl0091--biotransformer-metabolite-activity-and-embedding-follow-up)
  and the [follow-up report](../analysis/2026-09-14-tdi-som-representation-screen.md#metabolite-activity-and-embedding-follow-up).

- **Four-endpoint regression placement reaches ST-RAE .475974, rank20:** the
  strict challenge-disjoint PXR Chemprop embedding has
  generic pIC50 signal, but a 10% addition is unstable for 2C9 and worsens 3A4.
  Validated hybrid files instead retain the full planned 1A2/D6 placement and
  move 2C9/3A4 by 25%, 50% or 100% toward archived source moments. Even full
  movement is only .180/.202 pIC50 mean absolute change for 2C9/3A4. With the
  user's high-risk, high-information preference, keep all-full as the exploratory
  choice. Candidate `f87e9914` / submission `8ef21e06` improves every endpoint;
  the largest movement is D6 1.353707→.565520 and rank83→20. Macro improves
  .719349→.475974 and rank69→20, while unchanged rank correlations isolate the
  gain as population placement rather than molecular ordering. See
  [RL0080](2026-09.md#rl0080--strict-pxr-chemprop-regression-transfer-and-hybrid-placement)
  and [the report](../analysis/2026-09-14-regression-pxr-hybrid-placement.md).
  Recorded decomposition attributes 80.96% of the improvement to D6 and 16.34%
  to1A2; changing 2C9/A4 adds another .006568 macro. Freeze the four placements
  as the new test-time baseline and prioritize a new D6 ranking signal.
  A corrected parent-CV transfer of SMARTCyp, ATTNSOM and CYP-MAP finds that
  SMARTCyp improves a fixed frozen-feature control from macro ST-RAE
  .617895 to .614356, but a 10% addition already worsens the latest development
  ensemble .586931 to .587254. ATTNSOM and CYP-MAP do not improve macro, and
  none supplies a material D6 ordering gain. Retain the features but do not add
  them to the current regression candidate. See
  [RL0082](2026-09.md#rl0082--som-representations-for-direct-inhibition-regression)
  and [the report](../analysis/2026-09-14-regression-som-representation-screen.md).

- **Challenge-disjoint PXR Chemprop is a strong independent A4 representation,
  but not a live-Caruana addition:** old PXR checkpoints overlap four blind-test
  structures, so a strict three-seed model was retrained after excluding all CYP
  connectivities. Morgan plus its 384-d embedding reaches target-prior MCC
  .373/.407 on UMAP/scaffold, and a tree/SVM/LGBM blend reaches .382/.426.
  Direct PXR response predictions are weak. Chemprop/live-Caruana cross-fit
  blending worsens; a lower-correlation classical PXR/Admetica 20% addition has
  +.0107 descriptive MCC but a bootstrap interval crossing zero. Keep both as
  saved alternatives; no submission. See
  [RL0078](2026-09.md#rl0078--challenge-disjoint-pxr-transfer-into-a4-tdi) and
  [the report](../analysis/2026-09-14-classification-pxr-transfer.md).

- **Live A4 Caruana is highly compressible:** the six effective members have
  mean blind-test pairwise Spearman .702, including .963 between the two Jeremy
  PCA256 readouts. OOF-only prediction distillation reduces the ensemble to the
  historical stack plus Jeremy PCA256 TabICL with blind Spearman .99472 and only
  seven top-249 swaps each way. Four members leave two swaps, and dropping the
  2% Morgan tail leaves none. A 15:1 addition of pure radius-0/1 unhashed counts
  gives a tiny +.00307 descriptive MCC change but lowers AUROC/AP and has an
  interval crossing zero; the ordinary-plus-unhashed composite consistently
  hurts. Keep the two-member recipe clean, or the exact three-type unhashed
  version only as an LB probe; no submission. See
  [RL0077](2026-09.md#rl0077--live-a4-caruana-member-correlation-and-compression)
  and [the report](../analysis/2026-09-13-classification-caruana-compression.md).

- **A4 block fusion helps CV but not the live ensemble; FDA TDI transfer is
  domain-shifted:** radius plus classification-Chemprop improves both internal
  splits, but hurts when added to live Caruana. The checksum-validated primary
  FDA/Instem source supplies 621 A4 TDI labels; overlapping challenge labels
  agree only64.7%. Independent transfer improves UMAP but worsens scaffold,
  similarity gating changes no calls, and alert-centered exact environments fail
  the dual readout. Keep the data, not the model; no submission. See
  [RL0076](2026-09.md#rl0076--block-fusion-primary-fda-tdi-source-and-local-alert-environments)
  and [the analysis](../analysis/2026-09-13-classification-dual-split-priorities.md).

- **Dual-split audit finds a robust radius/unhashed representation gain:** blind
  D6 test support is higher than even the training leave-one-out ceiling, so no
  honest all-row split can exactly match it. Keep UMAP as the OOD stress test and
  use the existing scaffold folds as the closest grouped secondary readout. PLS
  and local residual correction are rejected. Endpoint-specific unhashed Morgan
  improves ordinary Morgan on both splits, but exceeds the live D6 ensemble by
  only .00503 with an interval crossing zero; no submission. See
  [RL0075](2026-09.md#rl0075--dual-split-small-data-and-unhashed-radius-follow-up)
  and [the analysis](../analysis/2026-09-13-classification-dual-split-priorities.md).

- **Support gating rejected; validation mismatch identified:** a fold-excluded
  router among ten saved D6 ranks scores projected MCC .08495 versus .09485 live.
  UMAP OOF is much harder than test: median nearest-train Tanimoto .299 versus
  .471, with93.4% versus42.0% below.4. Keep UMAP as the OOD stress test; add a
  prespecified grouped secondary split matched to the test support histogram for
  likely-LB ordering. See
  [RL0074](2026-09.md#rl0074--fold-excluded-d6-support-routing) and
  [the analysis](../analysis/2026-09-13-classification-d6-support-gate.md).

- **Hashed FP positions hide D6 test novelty:** test density is46.23 bits/molecule
  versus48.04 D6 and49.80 A4. Nearly all test bit positions occur in training, but
  unhashed IDs show4.58% occurrence-weighted novelty against D6 and3.03% against
  A4. D6 median nearest-train Tanimoto is only.471 and42% of test is below.4.
  Treat this as structural coverage/frequency shift, not a fingerprint alignment
  bug; favor larger/unhashed fingerprints and support-aware blends. See
  [RL0073](2026-09.md#rl0073--traintest-morgan-density-and-environment-coverage)
  and [the analysis](../analysis/2026-09-13-classification-fp-train-test-shift.md).

- **Morgan attribution separates the endpoint regimes:** D6 has only two
  non-ubiquitous FDR-significant bits and gains from radius-2-only environments,
  while A4 has 110 significant bits and gains from radius0/1. D6's top-50
  coefficients alone fail and their removal improves AUROC, so do not turn hashed
  top bits into manual alerts. Next use prespecified radius-isolated and unhashed
  collision-control features; no live change. See
  [RL0072](2026-09.md#rl0072--morgan-bit-attribution-and-radius-ablation) and
  [the analysis](../analysis/2026-09-13-classification-fp-attribution.md).

- **Fixed-feature model variation finds a small D6 complement:** with Morgan
  radius-2 bits and 49 calls fixed, honest fold-excluded selection chooses RBF-SVM
  or ComplementNB and moves projected MCC from .09485 to .09988. The +.00503
  conditional interval crosses zero; the validated candidate swaps only six live
  calls each way and was not uploaded. Keep the current live reference and retain
  this only as a low-risk optional LB probe. See
  [RL0071](2026-09.md#rl0071--fixed-fingerprint-classical-d6-model-comparison) and
  [the report](../analysis/2026-09-13-classification-d6-fixed-fp-models.md).

- **Broader classical D6 screen is promising but unstable:** radius-3 Morgan is
  the best pooled extension and FP+Chemprop wins three of five fold-specific
  selections, but fold-excluded model choice does not beat the live D6 recipe.
  Keep the current submission; next test only a prespecified radius2/radius3/torsion
  late-rank ensemble. See
  [RL0070](2026-09.md#rl0070--classical-fingerprint-and-chemprop-d6-screen) and
  [the report](../analysis/2026-09-13-classification-d6-fp-chemprop.md).

- **Diverse classification candidate reaches macro MCC0.352961, rank13:** the
  user-approved logistic-blend49 D6 plus Caruana-rank249 A4 improves D6
  .176698→.270958 and A4.375036→.434963. Macro rises from.275867 at rank48 and
  exceeds the public TabICL baseline.324967. Gain decomposition attributes61.1%
  to the D6 blend,36.1% to restoring the A4 call rate, and2.8% to A4 Caruana
  reranking. See
  [RL0069](2026-09.md#rl0069--classification-lb-gain-decomposition) and
  [the result analysis](../analysis/2026-09-13-classification-diverse-lb-result.md).

- **A materially new classification candidate is prepared:** a2:1 rank blend of
  currentD6 and weak Morgan logistic regression uses49 D6 calls, removing14 and
  adding20 relative to the live vector. Projected hidden-prior MCC improves on all
  five folds; the post-screen conditional interval is narrowly positive. Combined
  with the249-call Caruana A4 rerank, both endpoints change. Candidate validates;
  no upload. See
  [RL0066](2026-09.md#rl0066--diverse-cyp2d6-rankers-after-lb-prior-analysis) and
  [the report](../analysis/2026-09-13-classification-d6-diverse-models.md).

- **Classification LB correction favors endpoint rollback, not a global rate
  shift:** inferred hidden positive rates are about7.06% forD6 and29.16% forA4.
  The current D6 contraction helped, while the A4 expansion hurt. CurrentD6 plus
  oldA4 has known macro MCC0.303676 versus0.275867 current, though still below the
  public TabICL baseline0.324967. Prior-weighted LightGBMs and D6 rank blends do
  not justify replacing currentD6; an A4 Caruana rerank is retained only as a
  higher-risk test. No upload. See
  [RL0065](2026-09.md#rl0065--classification-leaderboard-prior-correction) and
  [the analysis](../analysis/2026-09-13-classification-lb-prior-correction.md).

- **Direct TDI structural-alert concatenation rejected:** a fixed 12-category
  bioactivation SMARTS panel contains real univariate signal, including 2.64-fold
  CYP2D6 methylenedioxyphenyl enrichment, but worsens honest inner-threshold macro
  MCC for both Morgan and Morgan+Mordred LightGBM. Retain only as an interpretable
  diversity feature for future hypotheses; do not adopt in the current candidate.
  See [RL0064](2026-09.md#rl0064--explicit-tdi-structural-alert-augmentation) and
  [the comparison](../analysis/2026-09-13-classification-tdi-structural-alerts.md).

- **Large-effect regression placement candidate prepared:** current live CYP2D6
  has ST-RAE1.353707, Spearman.431249 and R²−.813077. An exact current-vector
  candidate shifts only CYP2D6 from mean/SD4.6559/.3993 to the leaderboard-informed
  3.57/.90, preserving all ranks and the other three endpoints. Macro below.60
  requires placed D6 below.8763; the public source's stronger vector scored.5653
  at this placement, but that result is not assigned to ours. The combined
  CYP1A2+CYP2D6 version is the designated next regression candidate and passes
  submission dry run. See
  [RL0063](2026-09.md#rl0063--large-effect-lb-informed-cyp2d6-placement-candidate).

- **Regression ensemble pruned target by target:** the preferred simplified
  result uses3/3/2/5 top-level components for1A2/2C9/2D6/3A4 and compresses the
  internal3A4 Caruana from20 LightGBMs to two. Caruana is removed from the other
  three endpoints, so its total falls77→2; cross-fitted macro ST-RAE is0.594440
  versus0.596414 for the submitted five-way
  mean. The candidate submitted on September 13 at 05:51 UTC scored leaderboard
  macro ST-RAE0.719349 (rank67), improving the prior0.721167 by0.001817. See
  [RL0061](2026-09.md#rl0061--cyp1a2-caruana-compression-and-removal)
  and [the pruning report](../analysis/2026-09-13-regression-component-pruning.md).

- **Classification Caruana and compact candidates completed:** 18 candidates,
  nine selection procedures; 10 new triple-excluded auxiliary models and10
  historical inner readouts. D6 best tested Caruana MCC0.084431 <fixedmean0.096228;
  A4 sparseCaruana0.356736 <historicalTabPFN0.365149. FixedMorgan D6single scores
  0.093546; exact A4 restoration and compact candidates prepared/dry-run validated,
  no upload. See [RL0057](2026-09.md#rl0057--classification-caruana-and-compact-candidates)
  and [completed comparison](../analysis/2026-09-12-classification-caruana.md).

- **Comprehensive D6 classification audit completed:** 70 new controlled fits,
  raw-call reconciliation, actual official scores, 25-model probability audits,
  structural coverage and assay-cohort controls. Fewer labels alone do not explain
  the gap. A4 includes1,250 incomplete-potency-pair negatives versus4 forD6;
  original UMAP is structurally harsher than test even with identical neighbor
  pools. No bulk label/identity failure found; no D6-specific noise ceiling
  established. See [RL0056](2026-09.md#rl0056--comprehensive-d6-classification-explanation-audit)
  and [full audit](../analysis/2026-09-12-d6-comprehensive-audit.md).

- **Classification result confirmed:** MA-MCC0.275867 (rank45), D6 MCC0.176698
  improves from0.001644, while A4 MCC0.375036 falls from0.430653. Preserve the
  old A4 baseline for a future endpoint combination. See
  [RL0055](2026-09.md#rl0055--classification-ensemble-leaderboard-result).

- **Classification submitted:** D6 broad mean and A4 balanced LGBM mean from the
  completed round, 43/378 test positives. Official intake confirmed at September
  12 12:27 UTC; scoring pending. See [RL0054](2026-09.md#rl0054--classification-round-submission).

- **Classification round completed during regression cooldown:** 17 individual
  recipes and seven fixed probability means on the original UMAP folds. New
  inner-threshold MCC reaches 0.096228 for D6 and 0.338112 for A4, from different
  means. Keep the historical A4 baseline, whose AUROC remains higher. Models,
  inner/outer predictions and test probabilities saved; no upload. This resumes
  classification after its earlier deferral. See [RL0053](2026-09.md#rl0053--classification-round-during-regression-cooldown)
  and [completed comparison](../analysis/2026-09-12-classification-round.md).

- **Calibration comparison complete:** 30 fixed methods; native macro improves
  0.622574 to 0.601682, but D6 only 0.823983 to 0.819776 with a conditional
  uncertainty interval crossing zero. Source placement is separately applied to
  the exact previously submitted vector; D6-only CSV changes no other endpoint.
  Files validated, no upload. See [RL0051](2026-09.md#rl0051--broad-calibration-screen-and-same-vector-placement-controls)
  and [completed comparison](../analysis/2026-09-12-placement-round.md).

- **Native community learning reproduction complete:** all 25 fits use source
  data and native Chemprop defaults. Fixed mean gives D6 0.823983 / Spearman
  0.472811, essentially tied with previous 0.825436 / 0.472752; macro worsens
  to 0.622574. A separate audit identifies different local ST-RAE denominators,
  not a demonstrated live-LB scoring bug. No adoption or submission. See
  [RL0050](2026-09.md#rl0050--native-community-learning-reproduction-and-metric-definition-audit)
  and [completed comparison](../analysis/2026-09-12-scp-training-reproduction.md).

- **2D6-focused follow-up complete:** specialist scope/seed means, external efficacy
  features, actual joint public graph learning, interval-width decisions and kernel
  regressors. Fixed scope-mean addition gives D6 0.822232 versus 0.825436, but the
  bootstrap interval includes zero; expanded Caruana is slightly worse (0.825952).
  No submission. See [RL0049](2026-09.md#rl0049--cyp2d6-specialists-joint-external-responses-and-interval-decisions)
  and [completed program](../analysis/2026-09-12-d6-focus-round.md).

- **Report-derived round complete:** TDI auxiliaries, Hill-linked response heads,
  geometry/charge proxies, external substrate scores and joint descriptors.
  Expanded 126-candidate weight-crossfit macro is 0.582816 versus 0.583257.
  Small 1A2 gain; 2D6 essentially unchanged. Atom-level SoM remains untested.
  No submission. See [RL0048](2026-09.md#rl0048--report-derived-auxiliary-and-descriptor-experiments)
  and [completed program](../analysis/2026-09-12-report-experiments.md).

- **Global Caruana completed:** 107 saved candidate identities across model families.
  Weight-crossfit macro is 0.583257 versus 0.586931; gains are in 2C9 and 3A4,
  while 1A2 and 2D6 slightly worsen. Full-OOF fitted score 0.573705 is a fitting
  result, not independent CV. Weights and predictions are retained; no submission.
  See [RL0047](2026-09.md#rl0047--global-caruana-across-saved-model-families)
  and [the global comparison](../analysis/2026-09-12-global-caruana.md).

- **Jazzy comparison complete:** standard ETDG features for all 6,896 compounds
  take about two minutes; 12 net-charged structures retain missing values. Paired
  LGBM/TabPFN controls show small mixed gains, with no improvement over the latest
  ensemble. Retain features and results; latest candidate and submission remain
  unchanged. See [RL0046](2026-09.md#rl0046--dg-jazzy-descriptors-and-controlled-readouts)
  and [the controlled comparison](../analysis/2026-09-12-jazzy-readout.md).

- **Latest focused round:** prioritize1A2/2D6, allowing useful auxiliary tasks.
  Improved log2fc features and a complementary2D6 readout yield OOF-only
  candidates around0.5869; no new submission from that round. The later simplified
  regression ensemble is the latest actual LB at0.719349.
  See [RL0045](2026-09.md#rl0045--focused1a22d6-experiments-and-improved-primary-features)
  and [the focused program](../analysis/2026-09-12-focus-program.md).

- **2026-09-12 update:** The regression five-component ensemble reached development
  OOF ST-RAE **0.596414** on the retained parent CV. Its full-data refit has now
  been submitted without calibration; its scored LB is 0.721167 (RL0044).
  See [RL0043](2026-09.md#rl0043--full-data-five-component-submission) and the
  [deployment report](../analysis/2026-09-12-full-ensemble-submission.md).
  This supersedes the older model/pause/split status bullets below, which retain
  links to the earlier decisions. Regression remains the priority.

- Interval-aware TabPFN/LightGBM and neighborhood trials are complete. The best
  numerical blend scores 0.688853, but its incremental gain over the simple q0.55
  ensemble is uncertain; retain alternatives without replacing the submitted baseline
  ([RL-0017](2026-09.md#rl-0017--interval-aware-readouts-and-neighborhood-alternatives)).
- Continued regression trials favor a simple three-member TabPFN 2.6 q0.55
  ensemble in exploratory OOF (ST-RAE 0.691191). Local q0.55/median candidates
  are saved, with hidden-test population shift unresolved
  ([RL-0016](2026-09.md#rl-0016--distribution-summaries-and-simple-ensembles-improve-regression-cv)).
- External activity modeling is paused. A PXR-inspired calibration screen supports
  modest interval-loss gains but does not identify the hidden test-population
  correction. Keep full/half-strength maps for prospective evaluation; preserve
  the existing submission
  ([RL-0015](2026-09.md#rl-0015--pxr-inspired-regression-calibration-without-external-activity-data)).
- Focus on regression; defer classification tuning until regression stalls or
  priorities change. Saved-prediction diagnosis is complete
  ([RL-0012](2026-09.md#rl-0012--regression-diagnosis-and-existing-chembl-inventory)).
  ChEMBL 36 and PubChem auxiliaries are now curated with source lineage and
  challenge-connectivity exclusion: 24,075 external compounds, 93,119 labels
  ([RL-0013](2026-09.md#rl-0013--external-chembl-and-pubchem-auxiliary-dataset)).
- The first matched neural comparison is complete: external auxiliaries improve
  the point estimate of MA-ST-RAE from 0.824695 to 0.797301, but do not reach the
  adopted TabPFN's 0.724704. Retain the dataset and baseline; do not adopt this
  neural recipe. The macro conditional interval slightly crosses zero
  ([RL-0014](2026-09.md#rl-0014--first-matched-external-data-regression-comparison)).
- Pooled OOF means are nearly calibrated, but CYP2D6 ranking remains weak across
  folds. Hidden-test calibration is unresolved; retain UMAP and separate future
  calibration experiments from the external-data comparison.

- Baseline-v1 uses the augmented TabPFN 2.6 fold ensemble for **both** tracks by
  user decision. First LB scores are recorded: regression ST-RAE 0.868726,
  classification MCC 0.216149; CYP2D6 is the main weakness
  ([RL-0011](2026-09.md#rl-0011--first-baseline-leaderboard-scores)).

- Use the existing UMAP five-fold assignment for routine experiments on both tracks.
  No routine random or scaffold evaluation. Reconsider the split only if observed
  leaderboard results diverge substantially from CV; record the evidence first.
- CheMeleon-initialized auxiliary pretraining is complete ([RL-0005](2026-09.md#rl-0005--chemeleon-initialized-log2fc-pretraining)).
  The first fold-matched embedding-only TabPFN readout is complete ([RL-0006](2026-09.md#rl-0006--embedding-only-tabpfn-readout)).
  Regression improves over the initial LightGBM pipelines; classification MCC does
  not. Predicted-log2fc readout inputs require nested cross-fitting. Scratch
  pretraining remains deferred.
- The PCA256 + predicted log2fc4 TabPFN experiment is complete ([RL-0007](2026-09.md#rl-0007--foundation-pca-and-predicted-log2fc-with-tabpfn)).
  Changes versus the embedding readout are marginal; CYP2D6 TDI remains weak.
- Mordred augmentation gives a small regression improvement; the first direct-FT
  recipe underperforms ([RL-0008](2026-09.md#rl-0008--mordred-augmentation-and-direct-chemeleon-ft)).
  Classification remains below the earlier LightGBM MCC.
- Retain TabPFN 2.6 for the augmented feature recipe: v3 with eight or nine
  estimators underperforms on both tracks ([RL-0009](2026-09.md#rl-0009--tabpfn-30-on-regression-and-classification)).
- See [RL-0003](2026-09.md#rl-0003--umap-only-and-local-research-log) for the current decision.

## Recording policy

Append dated entries to `YYYY-MM.md`, with stable sequential IDs (`RL-0001`, etc.).
Use the actual event date/time with an explicit timezone when available. Label
retrospective entries as such; do not invent timestamps. One entry may cover a
coherent batch of runs. Include scientifically informative negative and discarded
results, not only improvements. Do not turn path mistakes, schema typos, preflight
failures, or empty interruptions into research findings. Substantive experiments go
through `pixi run cyp run`; small read-only inspections and preflight checks do not
require a run. Failed runs with no usable scientific result may be removed after a
concise reason is recorded. Link the log entry in a retained run note when useful.

Keep each entry short and evidence-backed:

- Question / hypothesis and track(s).
- Method / changes, configuration and all relevant run IDs.
- Results and limitations; link full reports or artifacts instead of copying tables.
- Interpretation, decision (keep/discard/inconclusive), and next step.

Never silently rewrite an earlier conclusion. Append a new entry referencing the
superseded one. Small factual corrections should be explicitly marked. Keep this
index's current direction up to date with links to the decisions supporting it.

PostgreSQL stores structured run, metric, split and prediction records. Run folders
store immutable source/config/input snapshots, models and outputs. This log stores
the reasoning connecting those records. Detailed reports live in `docs/analysis/`;
official source text remains separate in `references/official/`.

GitHub Issues are optional for bounded implementation tasks or actionable bugs;
they are not a second research log. Existing `--issue` run metadata stays optional.
Reports and research entries default to English; conversation language is independent.

## Entries

| Entry | Topic |
|---|---|
| [RL-0001](2026-09.md#rl-0001--preprocessing-and-descriptor-exploration) | Retrospective: preprocessing and descriptor exploration |
| [RL-0002](2026-09.md#rl-0002--initial-split-baseline-comparison) | Retrospective: LightGBM split comparison |
| [RL-0003](2026-09.md#rl-0003--umap-only-and-local-research-log) | UMAP-only evaluation and local research log |
| [RL-0004](2026-09.md#rl-0004--frozen-chemeleon-raw-2048-fingerprints) | Frozen CheMeleon features and checkpoint provenance |
| [RL-0005](2026-09.md#rl-0005--chemeleon-initialized-log2fc-pretraining) | Five-fold auxiliary pretraining from CheMeleon |
| [RL-0006](2026-09.md#rl-0006--embedding-only-tabpfn-readout) | First frozen embedding TabPFN readout |
| [RL-0007](2026-09.md#rl-0007--foundation-pca-and-predicted-log2fc-with-tabpfn) | PCA256 + nested predicted log2fc4 with TabPFN |
| [RL-0008](2026-09.md#rl-0008--mordred-augmentation-and-direct-chemeleon-ft) | Mordred-only, augmentation and direct FT |
| [RL-0009](2026-09.md#rl-0009--tabpfn-30-on-regression-and-classification) | TabPFN 3.0 with eight and nine estimators |
| [RL-0010](2026-09.md#rl-0010--tabpfn-baseline-v1-preparation-and-submission) | Baseline v1, private-profile upload automation and receipts |
| [RL-0011](2026-09.md#rl-0011--first-baseline-leaderboard-scores) | First baseline LB results, endpoint breakdown and CV comparison |
| [RL-0012](2026-09.md#rl-0012--regression-diagnosis-and-existing-chembl-inventory) | Regression diagnosis, classification deferral and local ChEMBL inventory |
| [RL-0013](2026-09.md#rl-0013--external-chembl-and-pubchem-auxiliary-dataset) | ChEMBL/PubChem acquisition, curation and source lineage |
| [RL-0014](2026-09.md#rl-0014--first-matched-external-data-regression-comparison) | Matched neural external-data comparison and conditional uncertainty |
| [RL-0015](2026-09.md#rl-0015--pxr-inspired-regression-calibration-without-external-activity-data) | PXR-inspired interval calibration, label-isolation checks and unresolved test shift |
| [RL-0016](2026-09.md#rl-0016--distribution-summaries-and-simple-ensembles-improve-regression-cv) | Distribution summaries, simple ensembles and local prospective candidates |
| [RL-0017](2026-09.md#rl-0017--interval-aware-readouts-and-neighborhood-alternatives) | 120 interval/point readout fits, neighborhood alternatives and small incremental gains |
| [RL0091](2026-09.md#rl0091--biotransformer-metabolite-activity-and-embedding-follow-up) | Fold-safe metabolite activity scoring and frozen parent/product embedding follow-up |
| [RL0092](2026-09.md#rl0092--classification-simplification-raw-embeddings-and-peer-replacement) | Classification hierarchy, raw/PCA comparison, compression and peer replacement |
| [RL0093](2026-09.md#rl0093--extensible-classification-stacking-and-safe-anchoring) | Extensible nested stacking, expert-set trial and safe anchoring |
| [RL0094](2026-09.md#rl0094--d6-rawpca-block-prioritization-and-remaining-a4-raw-controls) | D6 raw/PCA candidate priority and untested A4 raw controls |
| [RL0095](2026-09.md#rl0095--a4-matched-raw-readouts-and-fixed-mean-simplification) | A4 raw600 controls and fixed-mean simplification |
| [RL0096](2026-09.md#rl0096--d6-mean-block-plus-a4-equal-seven-leaderboard-result) | D6 raw/PCA block plus A4 equal-seven leaderboard result |
| [RL0097](2026-09.md#rl0097--regression-chemeleon-raw-readout-and-full-member-placement) | Regression CheMeleon raw2048 readout and full-member placement |
| [RL0098](2026-09.md#rl0098--project-owned-external-encoder-and-som-representation-screen-for-d6) | Project-owned external encoder and SoM representation screen for D6 |
| [RL0099](2026-09.md#rl0099--corrected-all-target-jeremy-replacement-evaluation) | Corrected all-target project-owned versus Jeremy encoder comparison |
| [RL0100](2026-09.md#rl0100--jeremy-medium-public-source-family-acquisition) | Public source-family acquisition for a project-owned Jeremy-medium replacement |
| [RL0101](2026-09.md#rl0101--project-owned-encoder-v2-closes-the-jeremy-medium-gap) | Project-owned encoder v2 replacement and complement screen |
| [RL0102](2026-09.md#rl0102--focused-intrinsic-cyp-task-objective) | Focused CYP continuous pretraining objective and intrinsic confirmation |
| [RL0103](2026-09.md#rl0103--focused-task-intrinsic-gains-do-not-transfer-to-raw-lightgbm) | Focused-task raw LightGBM transfer check and negative result |
| [RL0104](2026-09.md#rl0104--downstream-guided-task-objectives-recover-owned-transfer) | Downstream-guided task-objective screen and held-fold confirmation |
| [RL0105](2026-09.md#rl0105--full-oof-task-design-audit-and-follow-ups) | Full-OOF task/source audit and target-family follow-ups |
| [RL0106](2026-09.md#rl0106--jeremy-medium-gap-points-to-joint-multitask-optimization) | Jeremy checkpoint, overlap and representation-geometry reanalysis |
| [RL0107](2026-09.md#rl0107--joint-sparse-matrix-training-nearly-closes-the-owned-gap) | Project-owned joint sparse-matrix training and target-specific transfer |
| [RL0108](2026-09.md#rl0108--small-joint-encoder-refinements-do-not-surpass-jeremy) | Target weighting, nested readout and weight-space interpolation follow-up |
| [RL0109](2026-09.md#rl0109--fold-safe-challenge-train-and-auxiliary-adaptation) | Fold-safe direct pIC50, primary log2FC and direct Emax adaptation |
| [RL0110](2026-09.md#rl0110--jeremy-representation-distillation-surpasses-the-frozen-teacher) | Jeremy raw600 representation distillation and external-task anchoring |
| [RL0111](2026-09.md#rl0111--jeremy-initialized-pubchem-only-adaptation-is-marginal-for-d6) | Jeremy-initialized PubChem-only CYP2D6 adaptation |
| [RL0112](2026-09.md#rl0112--predictor-views-and-flat-checkpoint-ensembles-improve-d6) | Predictor views, diverse readouts and target-specific flat checkpoint banks |
| [RL0113](2026-09.md#rl0113--three-full-fit-d6-checkpoint-submission-candidates) | Full-fit stable-five, equal-ten and checkpoint-only D6 candidates |
| [RL0114](2026-09.md#rl0114--tabpfn-35-target-specific-flat-banks-and-all-target-candidates) | TabPFN 3.5 target-specific flat banks and all-target submission candidates |
| [RL0115](2026-09.md#rl0115--continuous-tdi-representations-improve-the-3a4-regression-bank) | Continuous TDI representations in the target-specific 3A4 regression bank |
| [RL0116](2026-09.md#rl0116--heterogeneous-lightgbm-feature-combinations-restore-tabular-diversity) | Heterogeneous LightGBM feature combinations and restored tabular diversity |
| [RL0117](2026-09.md#rl0117--correlation-audited-regression-banks-remove-near-duplicate-peers) | Correlation-audited regression banks and duplicate removal |
| [RL0118](2026-09.md#rl0118--regression-adopts-a-parent-retained-dual-split-readout) | Parent-retained dual-split regression readout |
| [RL0119](2026-09.md#rl0119--foundation-readouts-and-lgbm-swaps-produce-v8-submission-candidate) | Foundation readouts, LGBM swaps and v8 candidate |
| [RL0120](2026-09.md#rl0120--submission-aware-placement-produces-an-adaptive-d6-v9-candidate) | Submission-aware adaptive-D6 placement v9 |
| [RL0121](2026-09.md#rl0121--adaptive-d6-v9-selected-as-the-next-regression-submission-candidate) | Adaptive-D6 v9 selected as next regression candidate |
| [RL0122](2026-09.md#rl0122--public-classification-reports-leave-one-bounded-d6-tree-ensemble-test) | Public classification report review and bounded D6 tree-ensemble follow-up |
| [RL0123](2026-09.md#rl0123--public-tree-recipe-and-regression-feature-d6-transfer-do-not-clear-the-operating-boundary) | D6 public tree reproduction, regression-feature transfer and fixed-rank audit |
| [RL0124](2026-09.md#rl0124--newly-released-tabpfn-35-improves-d6-ranking-but-ties-the-operating-boundary) | TabPFN 3.5-fast/full D6 classification screen and blend audit |
| [RL0125](2026-09.md#rl0125--tabpfn-35-also-improves-a4-cv-but-does-not-displace-the-rank-1-vector) | TabPFN 3.5-fast/full A4 symmetry check and equal-seven blend audit |
| [RL0126](2026-09.md#rl0126--superseded-deferred-a4-tabpfn-35-manual-blend) | Superseded deferred A4 TabPFN 3.5 manual blend |
| [RL0127](2026-09.md#rl0127--restore-the-no-manual-blend-rule-for-new-classification-experts) | Withdraw manual TabPFN blend candidates and restore common-bank selection |
| [RL0128](2026-09.md#rl0128--bank-level-tabpfn-35-addition-and-estimator-stabilization-both-fail-d6) | D6 TabPFN 3.5 bank extension and 32-estimator stability rejection |
| [RL0129](2026-09.md#rl0129--test-like-validation-reverses-the-d6-tabpfn-conclusion) | Random/test-support D6 validation and equal-vote TabPFN probes |
| [RL0130](2026-09.md#rl0130--support-selected-stratified-random-folds-expose-large-d6-split-variance) | Support-selected stratified-random D6 folds and multiseed TabPFN audit |
| [RL0131](2026-09.md#rl0131--support6-d6-consensus-submitted-as-a-prospective-leaderboard-probe) | Submitted support6 D6 consensus leaderboard probe |
| [RL0132](2026-09.md#rl0132--support6-d6-leaderboard-probe-fails-decisively) | Support6 D6 leaderboard failure and candidate rejection |
| [RL0133](2026-09.md#rl0133--label-aware-neighborhoods-explain-the-support6-d6-failure) | Label-aware neighborhood and correlated-vote diagnosis of the D6 failure |
| [RL0134](2026-09.md#rl0134--fail-closed-classification-candidate-preflight) | Fail-closed classification candidate safety gate |
| [RL0135](2026-09.md#rl0135--regression-style-heterogeneous-feature-combinations-for-d6) | Expanded D6 heterogeneous feature combinations and controlled add/swap audit |
| [RL0136](2026-09.md#rl0136--regression-trained-jeremy-anchoring-transfers-to-d6-classification) | Jeremy-anchored external adaptation with a new CheMeleon-backed D6 family |
| [RL0137](2026-09.md#rl0137--direct-d6-jeremy-anchoring-does-not-improve-the-representation) | Negative fold-safe D6-label Jeremy anchoring experiment |
| [RL0138](2026-09.md#rl0138--submission-aware-v9-regression-upload-and-its-leaderboard-verdict) | Submission-aware v9 regression upload, 1A2 gain and D6 spread rejection |
| [RL0139](2026-09.md#rl0139--new-anchored-jeremychemeleon-d6-bank-is-stronger-but-not-selection-stable) | New-family-only anchored Jeremy/CheMeleon D6 bank and scaffold confirmation |
| [RL0140](2026-09.md#rl0140--combined-oldnew-d6-selection-retains-the-new-family-without-legacy-support) | Combined legacy/new D6 pool, flexible anchor-support selection and legacy rejection |
| [RL0141](2026-09.md#rl0141--one-flat-oldjeremy-chemeleon-pool-does-not-yield-an-mcc-caruana-gain) | One flat 43-member D6 pool and held-fold Caruana comparison |
| [RL0142](2026-09.md#rl0142--weight-guided-pool-adjustment-improves-d6-without-a-member-cap) | Weight-guided pool adjustment and unrestricted D6 Caruana improvement |
| [RL0143](2026-09.md#rl0143--ap-mass90-full-fit-is-diverse-but-fails-enforced-preflight) | AP-mass90 full fit, seven-member weights and enforced preflight failure |
| [RL0144](2026-09.md#rl0144--blind-stable-ap-mass90-d6-challenge-candidate-submitted) | Blind-stable AP-mass90 D6 challenge submission with A4 unchanged |
| [RL0145](2026-09.md#rl0145--ap-mass90-leaderboard-result-confirms-d6-oof-to-blind-mismatch) | AP-mass90 D6 leaderboard loss and OOF transfer failure |
| [RL0146](2026-09.md#rl0146--rebuild-classification-evaluation-around-model-and-call-rate-selection) | Classification model and call-rate evaluation protocol |
| [RL0147](2026-09.md#rl0147--three-accepted-d6-submissions-expose-two-distinct-oof-transfer-failures) | Three exact D6 submission comparison and primary proxy audit |
| [RL0148](2026-09.md#rl0148--historical-lb-concordance-is-the-criterion-for-a-replacement-oof) | Historical LB concordance as the OOF design target |
| [RL0149](2026-09.md#rl0149--exact-historical-classification-benchmark-frozen) | Seven exact classification submissions and endpoint LB transitions |
| [RL0150](2026-09.md#rl0150--full-molecule-d6-oof-changes-decisions-but-split-instability-remains) | Complete-pool D6 OOF refits and split-sensitive decision analysis |
| [RL0151](2026-09.md#rl0151--seven-member-ap-replay-removes-the-apparent-d6-oof-gain) | Seven-member AP-mass90 complete-pool replay and decision correction |
| [RL0152](2026-09.md#rl0152--direct-full-pool-control-replay-still-misses-the-lb-reversal) | Direct AP versus September 15 full-pool replay and unresolved LB reversal |
| [RL0153](2026-09.md#rl0153--scaffold-refit-and-blind-domain-weighting-do-not-restore-lb-concordance) | Scaffold refits and rejected train-to-blind density weighting |
| [RL0154](2026-09.md#rl0154--extratrees-and-assay-route-sensitivity-remain-split-dependent) | ExtraTrees direct comparison and rejected assay-route scoring-mask hypothesis |
| [RL0155](2026-09.md#rl0155--conservative-controlap-rank-blend-remains-uncertain) | Conservative control/AP rank blend and uncertain full-pool decision |
| [RL0156](2026-09.md#rl0156--analog-expansion-selection-exposes-the-validation-population-gap) | Blind analog-expansion selection and sparse D6 validation coverage |
| [RL0157](2026-09.md#rl0157--forty-nine-d6-calls-are-a-historical-comparator-not-a-constraint) | D6 model/rate audit rejects treating 49 calls as fixed |
| [RL0158](2026-09.md#rl0158--select-the-september-15-classification-vector-for-the-interim-submission) | September 15 classification vector selected for the pre-interim submission after cooldown |
| [RL0159](2026-09.md#rl0159--endpoint-wise-leaderboard-best-regression-submission) | Endpoint-wise leaderboard-best regression submission and exact scored confirmation |
| [RL0160](2026-09.md#rl0160--resubmit-the-exact-september-15-classification-vector) | Exact September 15 classification vector resent and scored at its prior MCC |
| [RL0161](2026-09.md#rl0161--intermediate-chemeleon-states-show-task-dependent-standalone-signal) | Intermediate CheMeleon depth screen across both tracks |
| [RL0162](2026-09.md#rl0162--fixed-intermediatefinal-chemeleon-blends-improve-simple-readouts) | Fixed intermediate/final CheMeleon OOF blends |
