# External assay representation program

The user authorized sustained experiments, parallel agents and a goal of matched
regression OOF macro ST-RAE below 0.60. The initial reference is the fixed four-way
mean at 0.658841, with 6,384 non-anchor observations on the existing parent-retained
folds. The goal is prospective; no improvement is claimed until measured.

## Design

1. Rebuild external observations with assay identity, species, endpoint type,
   relations/censoring, hitcalls, quality flags and original provenance retained.
   Start from local ChEMBL36 and PubChem1851. Acquire independent Tox21 P450 and
   ToxCast NVS_ADME measurements where available; deduplicate original experiments.
2. Keep architecture width600/depth5 and split/seed fixed for initial comparisons:
   assay-specific continuous regression, within-assay ranking, binary activity,
   and mixed objectives. Binary auxiliary learning serves the regression goal;
   competition TDI work remains deferred.
3. Evaluate external frozen embeddings with matched downstream models, initially
   LightGBM and TabPFN2.6. Compare raw embeddings with training-fold-only PCA200.
   Next evaluate challenge adaptation on training-fold pIC50/log2fc and external
   tasks, separately from the frozen representation comparison.
4. Retain overlap measurements separately for assay agreement and later
   training-fold-only calibration. Initial external encoders exclude all challenge
   connectivity blocks. Never use outer validation challenge labels to calibrate
   external labels or select training checkpoints.
5. Retain failed trials, modest gains and negative results. Compare per-endpoint,
   pooled OOF and neighbor scores, fixed blends and prediction/residual correlations.
   Future ensemble weight fitting requires separate selection/evaluation handling;
   fitting weights on the scored OOF is not an unbiased validation result.

Source guidance: local CYP Deep Research report18 (data scarcity), prior PXR
Deep Research regression report and QSAR literature, and prior implemented
`docs/track1_explain/phase2_pairwise_assay_rank.md`. These motivate experiments;
reported literature gains are not assumed to transfer to CYP.

The earlier PubChem600-only model remains a narrow comparator, not a recreation
of Jeremy's broader public CYP/ADME `chemprop_medium` pretraining corpus.
