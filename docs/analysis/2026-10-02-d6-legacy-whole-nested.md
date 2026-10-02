# Complete historical ensemble under the new D6 evaluation procedure

The user asks whether the complete existing classifier has been evaluated under
the procedure used for the new CheMeleon candidate. It has not: three historical
Morgan atoms were refitted, but cached whole-ensemble scores were only retrospective
comparators. This continuation refits the **complete eleven-atom recipe**.
The 175-call candidate is preserved and unsent. Its earlier recommendation was
made without this essential whole-model comparison and is provisional pending
the complete results.

## Frozen comparison and its interpretation

Run `20261002T095304331805Z-fad722ec` freezes the September 15 model family and
composition before the new comparison results. Inputs include the same canonical
molecules, labels, original outer assignments and challenge-independent external
features/encoder. No CYP-trained weights, old PCA/scalers or old scores enter
new fitting.

The eleven atoms are five fixed-hyperparameter LightGBM classifiers (Morgan,
Mordred, their concatenation, external class embeddings, external mixed
embeddings), two TabPFN 2.6 views (external mixed PCA256/raw600), a TDI-only
two-head fine-tuned graph model, and Morgan logistic C=.05 / RBF SVM C=1 /
ComplementNB alpha=1. Their historical raw-core/ranked-tail composition is fixed.
The two TabPFN predictions average into one core member, so the seven-member
core contains five LightGBMs, that mean, and the graph model.

Both original UMAP and scaffold outer assignments are used, with seeds
20261002/20261017. Each outer-training pool splits into five bag models; each
bag fits three deeper cross-fits for its bootstrap-mean-MCC threshold. Each
current graph fit additionally selects epochs on a group holdout **inside its
own training set**, then resets and refits within that same training boundary.
Graph training preserves max 80 epochs, patience 10, encoder/head AdamW rates
1e-5/1e-4, weight decay 1e-4, dropout .1, two observed-class BCE heads and batch
size 64. TabPFN scalers/PCA fit only current observed-D6 training rows.
The complete legacy recipe therefore needs 400 ensemble fits (4 per bag x
5 bags x 5 outer folds x 2 seeds x 2 strategies), not frozen OOF regrouping.

After measuring runtime and resource availability, computation is partitioned
by the four predetermined split/seed blocks. The initial run supplies its first
complete UMAP/20261002 block; the other blocks run concurrently as
`20261002T101409635209Z-967b7499` (UMAP/20261017),
`20261002T101409621831Z-5bf54a17` (scaffold/20261002), and
`20261002T101409637266Z-41e6beb1` (scaffold/20261017).
The initial process will be interrupted after that complete block to avoid
duplicating the other three blocks. Its immutable scientific results remain
retained. The merge requires identical feature/input hashes and common model
settings, exactly 100 fits and five complete outer folds in each block, and
exactly 400 total fits. No block is chosen by its score. This is a computational
change, not a scientific protocol change.

The final decision averages the ten outer scores minus their independently
trained cuts, exactly as in the new candidate. No 49-call budget is imposed.
Full-pool rank transforms see the prediction batch's molecules but never their
class labels, matching the historical transductive score composition. The review
requires exact molecule, endpoint, fold, group and label agreement with the
new-candidate outer records before paired comparisons.

This compares **model recipes under a shared decision procedure**. It does not
recreate every old submission's original operating-point procedure. Historical
recipe discovery itself happened on development labels; its discovery cannot
retrospectively be made nested by refitting a now-frozen recipe. The new recipe
was likewise selected after development inspections. These limitations remain
explicit, rather than calling either recipe a newly independent test.

Three conclusions must remain separate:

1. Matched whole-recipe comparison can establish an observed-label improvement
   relative to the complete existing family, rather than just a simple learner.
2. Reproducing the direction of past leaderboard changes would support the
   evaluation procedure's relevance to deployment. A same-method model win alone
   does not establish this historical concordance.
3. Neither comparison observes the 2,823 missing weak-primary D6 class outcomes;
   neither can by itself resolve the population identifiability problem.

Composition audit `20261002T095502819966Z-d8c4fea6` reconstructs all 6,146 archived
UMAP predictions from the eleven source atom records with **zero error** against
the historical raw-core composition. This checks the formula used in the active
refits; it is not additional fitting or generalization evidence.

The official [TDI starter notebook](https://github.com/OpenADMET/CYP-Challenge-Tutorial/blob/main/notebooks/TDI_prediction.ipynb)
currently illustrates descriptors with LightGBM. This is the tutorial recipe,
not verification of the architecture of any official leaderboard entry. The
new CheMeleon readout is still a single molecular representation family even
though its two-seed five-fold policy has ten fitted models. Architecture alone
is not evidence of improvement over the stronger historical composite.

The reused preprocessing script's generic `representation.json` initially marked
`no_pretrained_assets=true` when no fresh CheMeleon `embedding_file` was configured.
That generic flag omits the legacy branch's explicitly declared pretrained
external encoders/features and TabPFN checkpoint. Source manifests/configuration
and the fitting code identify these assets; the merged record corrects the flag
without altering the source records or any fits. Future preprocessing metadata
also recognizes the legacy branch.

## Results

The complete refit is in progress. No interim fold is used to select a recipe,
change the threshold objective, choose a count, or justify upload. The recorded
review will compare complete equal-seed decisions, paired conditional group
intervals and diagnostic matched rank fractions. Diagnostics at 49/750 and
175/750 are comparisons only, not prescribed new submission counts.
