# Classical fingerprint and Chemprop CYP2D6 screen

## Conclusion

Alternative fingerprints and frozen Chemprop embeddings contain complementary D6
signal, but this broad screen does not identify a stable replacement for the live
logistic-blend endpoint. Morgan radius 3 is the best pooled extension at 49 calls.
FP+Chemprop models win several fold-excluded selections, but the preferred feature
combination changes by fold. Retain the predictions; do not submit this winner.

## Design

Run `20260913T122320230702Z-d08c8ee5` trains 61 low-capacity base models on the
fixed five UMAP folds and evaluates 123 base/current-blend ranks. Features cover
Morgan radius 1/2/3 bits and counts, chiral Morgan, RDKit paths, atom-pair and
torsion counts, MACCS, two frozen Chemprop embeddings, and early FP+Chemprop
concatenation. Logistic regression uses three regularization strengths.

Each outer fold is excluded from model fitting. Candidate selection for a heldout
fold additionally uses only the other four OOF folds, making it more conservative
than choosing the best pooled OOF row.

## Results

At 49 calls, the best pooled row is the 2:1 original-ensemble/Morgan-radius-3
logistic blend (`C=0.05`), with projected target-prior MCC 0.104965. The live
radius-2 recipe scores 0.094847 in the same local comparison. The +0.0101 result is
promising but selected from 123 ranks.

Fold-excluded top-1, top-3 and top-10 selection scores 0.06109, 0.08008 and
0.08008, all below the live recipe's pooled 0.09485. The top-1 choice changes
among radius-3-count+mixed-Chemprop, torsion+mixed-Chemprop, radius-3 bits, and
RDKit-path-7+mixed-Chemprop. Chemprop fusion wins three of five fold selections,
but no single representation transfers consistently. Full-data top-three selection
drops Chemprop in favor of radius-3 bits and torsion counts.

## Interpretation and decision

The small-data hypothesis is supported: low-capacity linear fingerprints remain
competitive diversity components. Radius 3 suggests useful longer structural
context, while count and torsion fingerprints recur. Chemprop helps conditionally,
but embedding-only models do not lead and broad early-fusion selection is unstable.

Keep the live D6 and A4 endpoints. The next defensible experiment is a small,
prespecified late-rank ensemble of radius-2, radius-3 and torsion linear scores,
tested once with nested selection rather than another broad screen.

