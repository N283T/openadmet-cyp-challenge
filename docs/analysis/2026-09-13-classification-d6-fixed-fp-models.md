# Fixed-fingerprint CYP2D6 model-family comparison

## Conclusion

With Morgan radius-2 2,048-bit features and 49 positive calls fixed, nonlinear
RBF-SVM and Complement Naive Bayes are the only model families that repeatedly add
value to the live D6 rank. The fold-selected point estimate improves projected
target MCC from 0.09485 to 0.09988, but its conditional interval crosses zero.
Retain the validated six-swap candidate; do not submit it as a confirmed upgrade.

## Design

Run `20260913T123045515718Z-0df22d54` holds representation, split and positive
count constant while comparing 47 standard classifiers:

- Logistic regression and linear SVM.
- RBF-SVM and precomputed Tanimoto-kernel SVM.
- Random Forest and ExtraTrees.
- Bernoulli and Complement Naive Bayes.
- Ridge, SGD and shrinkage LDA.

It evaluates 143 standalone, original-ensemble blend and live-D6 blend ranks. For
each heldout UMAP fold, model selection uses only the other four OOF folds.

## Results

Pooled OOF favors a 2:1 blend of live D6 with ComplementNB (`alpha=1`), projected
MCC 0.11011 versus 0.09485 live. In fold-excluded top-1 selection, RBF-SVM (`C=1`)
is chosen for two folds and ComplementNB for three. The cross-fitted result is
0.09988. Top-3 scores 0.09485; top-5 again scores 0.09988. Tanimoto SVM appears in
larger ensembles but is not the leading incremental component. Forest methods
occasionally improve AUROC or AP without improving the sparse 49-call boundary.

Audit/package run `20260913T123420289849Z-b75917f6` gives a 10,000 UMAP-group
conditional delta of +0.00503 with interval [-0.01905, 0.03087]. Three folds have
identical fixed-rate calls, one worsens and one improves. The deployment choice is
live D6 blended 2:1 with ComplementNB. It retains 43 of 49 live test calls and
swaps six in each direction. The classification file validates and keeps A4
exactly unchanged; no upload was made.

## Interpretation and decision

The model-side hypothesis is useful: binary fingerprints support generative-count
and nonlinear-margin views that differ from logistic regression. However, the
remaining incremental signal is small compared with the already successful live
D6 revision. The six-swap candidate is a reasonable low-risk future probe, not a
demonstrated improvement. Further tuning of dozens of classifiers on these same
folds would mainly increase selection bias.
