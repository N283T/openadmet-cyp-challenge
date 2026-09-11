# Admetica embeddings alone into TabPFN 2.6

User requested the embedding-only model, without additional feature columns.
Fit run `20260911T150557650787Z-4a460d9e` uses the matched CYP's frozen Admetica
encoder, mean-pools its 300-dimensional message-passing output before BatchNorm,
and fits TabPFN 2.6 on pIC50. No probabilities, logits, descriptors, predicted
log2fc or external PCA are added. TabPFN's own preprocessing remains enabled.
The encoder is not fine-tuned.

All four sets of 6,896 keyed molecular embeddings and 20 fitted TabPFN states
are saved. Each endpoint uses current parent-retained five-fold assignments,
eight TabPFN estimators, and seed 42 + fold. Audit
`20260911T150745346553Z-6a171ac2` checks exact identity, fold, target and confidence
bound agreement for 6,384 OOF observations and registers scores in PostgreSQL.
No test pIC50 prediction or submission is produced.

## Results

| Model / fixed equal blend | All OOF macro ST-RAE | Neighbor macro ST-RAE |
|---|---:|---:|
| Admetica emb only | .743910 | .819000 |
| Vanilla CheMeleon emb only | .721992 | .895810 |
| FP+Mordred only | .728820 | .858036 |
| Predicted log2fc only | .698059 | .940962 |
| FP+Mordred + predicted log2fc (previous pair) | .673634 | .847109 |
| Admetica + predicted log2fc | .668797 | .829057 |
| Admetica + FP+Mordred | .690170 | **.795527** |
| Admetica + vanilla CheMeleon | .695482 | .826051 |
| FP+Mordred + predicted log2fc + Admetica | **.661279** | .803584 |
| Current full TabPFN | .688696 | .854329 |
| Current full TabICL | .695341 | .834264 |

Equal means are formed from complete component predictions, not feature
concatenation. Weights are not optimized. The triple changes the prior pair's
weights from 1/2 each to 1/3 each plus 1/3 Admetica.

| CYP | Admetica alone | Vanilla emb alone | Previous pair | New triple |
|---|---:|---:|---:|---:|
| 1A2 | .800361 | .785295 | .757879 | .740490 |
| 2C9 | .682876 | .640578 | .561990 | .554762 |
| 2D6 | .938661 | .907882 | .877768 | .873043 |
| 3A4 | .553743 | .554213 | .496899 | .476822 |

The new triple improves all four endpoint ST-RAEs over the prior pair. Macro
MAE improves .533516 -> .526064, and macro Spearman .591016 -> .607562.
Admetica-alone pIC50 Spearman is .500500/.580233/.403476/.746693, compared with
unchanged classifier probabilities .422950/.437027/.332296/.593393 on the same
validation population.

## Prediction complementarity

Pearson correlation of Admetica embedding/TabPFN predictions:

| CYP | vs FP+Mordred | vs vanilla emb | vs predicted log2fc |
|---|---:|---:|---:|
| 1A2 | .743544 | .780810 | .725615 |
| 2C9 | .729369 | .789638 | .741644 |
| 2D6 | .704010 | .778545 | .703432 |
| 3A4 | .854398 | .891165 | .841138 |

Correlation increases relative to unchanged classifier outputs after fitting a
pIC50 readout. The benefit is established by these fixed blends, not by low
correlation alone. The encoder-only model is not the strongest standalone model,
but is retained as an ensemble candidate. Caruana remains deferred until a larger
candidate pool, per user preference.

The external checkpoint's original training-set overlap has not been audited in
this experiment. Results describe the current external-pretraining protocol and
are not a new blinded leaderboard result. Current submission/model pointers are
unchanged. Detailed per-fold, neighbor, MAE, rank and residual-correlation tables
are in the audit run's outputs.
