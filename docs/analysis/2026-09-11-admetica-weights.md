# Admetica CYP weights and the CheMeleon redundancy claim

User clarified that Caruana selection is already planned, deferred until the
candidate pool is larger. Do not prioritize reintroducing it as a new idea.

## What our correlation result establishes

Frozen vanilla CheMeleon -> TabPFN predictions correlate .801–.916 with
FP+Mordred -> LightGBM on current folds. This is prediction similarity between
two complete pipelines, not equality of embeddings and descriptors or proof
that CheMeleon is a bad fine-tuning initialization. In fact, their equal mean
improves macro ST-RAE from .728820 (FP+Mordred alone) to .704776. CheMeleon adds
some useful information. The three-way mean .674589 does not improve on the
FP+Mordred plus pred_log2fc pair .673634, but does not establish zero marginal
information under every weighting or training recipe.

Jeremy's report describes redundancy in that author's own ensemble. Treat this
as an empirical observation, not a general fact about CheMeleon pretraining.
Comparing initializations needs matched downstream training, not just the
frozen-embedding correlation above.

## Acquired repository and relevant models

User-authorized ghq checkout:
`/home/nagaet/ghq/github.com/datagrok-ai/admetica`

Commit: `64eacdb7aab5a91c83b6f114379455deb3744ee5`.
Metadata inspection run: `20260911T144657359857Z-32d7892a`.
It snapshots four standard CYP inhibitor checkpoints plus the separate challenge
checkpoint and relevant source files. SHA-256 values and inspected pickle
metadata are retained in `checkpoint-inspection.json`.

The standard inhibitor models in `ADMET/metabolism/` cover all four scored CYPs
(and CYP2C19 separately). Metabolism documentation attributes them to PubChem
AID1851 binary inhibitor/non-inhibitor data, approximately 13k compounds per CYP.
They are not pIC50 regression heads; the transferable component is the encoder.

All four inspected standard models contain Chemprop BondMessagePassing,
mean aggregation, hidden width 300, depth 3 and a BinaryClassificationFFN.
This differs from vanilla CheMeleon's width 2048/depth6. Actual checkpoints
contain Chemprop 2-style objects despite the root requirements pinning 1.6.1.
Using the repository-wide old environment verbatim is therefore unnecessary.

## Concrete compatibility verification

Run `20260911T144903295084Z-ac8cd6ac` successfully loads the four standard
checkpoints, extracts encoder state dictionaries with strict key/shape matching,
and computes finite 300-dimensional embeddings for ethanol, benzene and caffeine.
Encoder-only `.pt` files are saved in that run's outputs. No installed package
or global environment was changed.

Loading retains `weights_only=True` with explicit known legacy aliases:
BinaryAccuracyMetric/BinaryAUROCMetric map to current metric classes, old NumPy
module names map to current equivalents, and metadata AttributeDict maps to a
plain dict. The classification predictor is not reused for regression.

Failed attempts remain recorded:
- `20260911T144803374670Z-52dffec2`: old NumPy global name was not mapped.
- `20260911T144835961838Z-fd18c358`: restricted loader rejected SETITEMS on
  Lightning AttributeDict; plain-dict mapping resolved it.

This is a load/forward smoke check, not regression performance evaluation or
proof of bitwise agreement with the original Admetica prediction environment.
Challenge/external pretraining overlap has not been audited in this inspection.

## Separate challenge-specific checkpoint

`examples/cyp_challenge/models/multitask_pac50.ckpt` is a different asset. Its
metadata has a 2048-wide/depth6 encoder and a regression head. The accompanying
training code consumes challenge labels and the production function trains on
all processed compounds. It must not be used as an external-only initialization
for our held-out challenge CV. This is not a restriction on the standard
AID1851 models above.

The same training script also explicitly supports CheMeleon warm starts, and its
docstring attributes improved task sharing to that initialization. Thus this
repository itself does not support dismissing CheMeleon as a starting point.

## Next controlled comparison

Use the standard AID1851-trained encoders as alternative initialization assets,
with fresh regression heads. First compare frozen per-CYP 300-dimensional
representations on the existing folds; then compare regression or log2fc
fine-tuning with the same validation protocol. Audit source-data overlap before
interpreting CV. Keep vanilla CheMeleon as a control; no current model is replaced
on the basis of prediction correlation alone.

Sources:
- https://github.com/datagrok-ai/admetica/tree/64eacdb7aab5a91c83b6f114379455deb3744ee5/ADMET/metabolism
- https://github.com/datagrok-ai/admetica/blob/64eacdb7aab5a91c83b6f114379455deb3744ee5/ADMET/metabolism/metabolism.md
- https://github.com/datagrok-ai/admetica/blob/64eacdb7aab5a91c83b6f114379455deb3744ee5/examples/cyp_challenge/scripts/train_cyp_multitask.py
