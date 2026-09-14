# PXR transfer features for CYP3A4 TDI classification

## Question and leakage policy

PXR activation induces CYP3A4 expression, but it is not the same mechanism as
the direct-versus-NADPH-preincubation TDI endpoint. This work therefore treats
prior PXR activity only as an auxiliary chemical phenotype or representation,
not as a signed proxy label for TDI.

Exact standardized-connectivity audit `d6dad5a9` found no current blind-test
overlap in the prior PXR pEC50 train or released-test sources. The old PXR
single-concentration source did overlap four blind-test connectivities, so its
frozen Chemprop checkpoints were not used. Failed precursor `79952629` records
the duplicate-connectivity join validation error fixed in the completed audit.

All new classical PXR models and the replacement Chemprop pretraining exclude
every one of the 6,889 current CYP connectivities before fitting.

## Classical PXR and predicted-phenotype blocks

Run `5177d5d9` trained Ridge, ExtraTrees, Random Forest and LightGBM PXR pEC50
models on 2,750 challenge-disjoint external structures. LightGBM was strongest
in five-fold scaffold CV (Spearman 0.720, RMSE 0.788). These scores were combined
with Morgan bits, four frozen Admetica inhibition probabilities and six external
CYP-substrate scores in runs `068c749c` and `c510a24e`.

PXR scores alone are weak for A4 TDI. The principal compact gain comes from
Admetica; PXR adds a smaller split-dependent contribution and substrate scores
usually dilute it.

| A4 model | UMAP target-prior MCC | Scaffold target-prior MCC |
|---|---:|---:|
| Morgan LightGBM | 0.1863 | 0.3268 |
| Morgan + PXR | 0.2036 | 0.3159 |
| Morgan + Admetica | 0.2302 | 0.3593 |
| Morgan + PXR + Admetica | 0.2396 | 0.3546 |
| Morgan + PXR + Admetica + substrate | 0.2412 | 0.3515 |

The PXR/Admetica model has OOF/test Spearman 0.579/0.744 to live Caruana. A
descriptive 20% rank addition raises target-prior MCC from 0.44394 to 0.45467
and changes 16 blind top-249 calls (`1954e23e`). Leave-one-fold-out selection
chooses the same model and weight on all five held-out folds (`d5069b46`), but
the 3,000-replicate UMAP-group interval for the +0.01074 change is
[-0.00363, 0.02371]. Retain only as an uncertain, low-correlation LB probe.

## Strict PXR Chemprop representation

Run `519bbb46` retrained the prior two-head log2fc Chemprop architecture at
8.25 and 33 uM on 6,616 challenge-disjoint PXR structures with seeds 42--44.
Scaffold-holdout Spearman ranges from 0.611 to 0.645 for the 8.25 uM head and
0.660 to 0.699 for the 33 uM head. It exports two prediction means, two seed
standard deviations, one embedding uncertainty summary and the mean 384-d
embedding for all current compounds.

The two predicted PXR response values are only weak A4 features. The embedding
is substantially stronger, showing that the gain is from external molecular
representation learning rather than a simple PXR-activity proxy.

| A4 model | UMAP target-prior MCC | Scaffold target-prior MCC |
|---|---:|---:|
| Chemprop scalar predictions + LR | 0.2224 | 0.2208 |
| Chemprop embedding + LightGBM | 0.3376 | 0.3963 |
| Morgan + Chemprop embedding + LightGBM | 0.3732 | 0.4071 |
| Morgan + PXR + Admetica + Chemprop embedding | 0.3639 | 0.4040 |

Run `8e8402cb` verifies the fixed embedding across model families. ExtraTrees is
best on UMAP (0.3623), LightGBM is best on scaffold (0.3963), and Random Forest
and RBF-SVM are also competitive. Thus the result is not unique to one readout.
Model-blend runs `3d84567e` and `cfe6cc15` find a stable three-way rank mean:
50% Morgan/Chemprop LightGBM, 25% embedding ExtraTrees and 25% embedding
RBF-SVM reaches 0.3824 on UMAP and 0.4255 on scaffold.

## Relationship to live Caruana and decision

Despite its strong standalone dual-split result, the strict Chemprop family is
already correlated with live Caruana (about 0.81 OOF and 0.86 blind test for the
Morgan/Chemprop model). Descriptive blends `e43a4239` give at most a small MCC
increase. Honest fold-excluded weight selection in `71b0e4fc` is unstable and
reduces target-prior MCC from 0.44394 to 0.43319. Do not add this family to live
Caruana.

Retain the strict PXR-Chemprop three-way model as a strong independent A4
alternative and representation bank. Retain the lower-correlation classical
PXR/Admetica 20% Caruana addition only as an optional LB probe whose bootstrap
interval crosses zero. Do not package or submit either result from this analysis.

## Nested peer-selection follow-up

Run `41ccf0ec` creates matched outer, rotated-inner and fold-test predictions for
five strict Chemprop readouts on both classification endpoints. Failed precursor
`986b33df` records a NumPy-integer/Polars indexing error. Run `2e6c1f38` then
adds these readouts as peers to the complete existing Caruana pool rather than
appending a fixed late blend.

Full-OOF AUROC selection assigns 15.15% total weight to strict Chemprop members,
so the selector detects their standalone signal. That does not survive honest
outer prediction at the deployed fixed call rate. Audit `9a7f115d` gives:

| Nested recipe | Original target-prior MCC | PXR-extended MCC |
|---|---:|---:|
| Caruana AUROC | 0.43473 | 0.42090 |
| Caruana Brier | 0.43166 | 0.40398 |
| Caruana MCC | 0.40552 | 0.40860 |
| Sparse-3 MCC | 0.42244 | 0.41475 |
| Bagged sparse-3 | 0.43319 | 0.41783 |

Only the unrestricted MCC objective gains, by 0.00308, while all more stable
recipes worsen. Blind ranks remain very close (Spearman 0.997--1.000, zero to 11
top-249 swaps each way). This confirms that peer replacement and calibration-
aware selection do not rescue the high-correlation Chemprop addition. The prior
decision stands: keep it separate, and if an LB diversity probe is desired, use
the lower-correlation 80/20 live-Caruana / classical-PXR-plus-Admetica rank mix.

## Candidate priority update

The 80% live-Caruana / 20% Morgan-PXR-Admetica A4 rank blend is now the first
candidate for the next classification probe. Since endpoint leaderboard scores
are separate, it should be paired with the previously validated six-swap
live-D6/ComplementNB rank blend rather than repeat the exact live D6 vector.

The strict PXR-Chemprop readouts were already eligible in the nested selection
bank for CYP2D6. Their addition does not improve the principal outer recipes:
MCC-Caruana changes from 0.08870 to 0.08486 and AUROC-Caruana from 0.09259 to
0.08432. Any D6 value from the PXR embedding would be generic representation
transfer rather than a direct PXR/CYP2D6 mechanism. It is therefore excluded
from the first candidate. No package or submission was created.
