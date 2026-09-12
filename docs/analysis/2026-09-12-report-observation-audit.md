# Independent audit of the observation-surrogate models

Model run: `20260912T093601580048Z-ec2b7210`.
Successful audit: `20260912T094310174909Z-8df31bae`.

All six methods contain exactly 6,384 canonical OOF rows. Compound, endpoint,
fold, split, labels, confidence bounds, and neighbor flags match the existing
`ae562a59` candidate. Predictions are finite. All 24 official endpoint scores
reproduce the source metrics within absolute tolerance 1e-12.

All 34,480 membership rows (6,896 compounds times five folds) match the canonical
parent assignments and label availability from `271e3020`. Outer training and
validation are disjoint. Inner training and validation partition outer training,
are disjoint by molecular group, retain anchors in training, and exactly reproduce
the seeded group split. Unassigned test compounds never enter training.

For all 15 checkpoints, every stored feature mean/standard deviation, head
mean/standard deviation, and direct-loss denominator exactly equals a fresh
float32 calculation using only the recorded outer training compounds. Model
normalization buffers match their saved arrays. The 647 input columns reproduce
the frozen 600-dimensional embedding plus 47 external scores in saved order.

All selected epoch counts reproduce the inner-validation history's improvement
criterion (1e-6 minimum improvement), and each outer refit has precisely that
many epochs with no outer-validation metric recorded. The selected counts are
1–5 epochs; this is an actual outcome of the early-selection protocol, not a
missing-refit or checkpoint-count error.

Independent CPU linear/ReLU reconstruction checks all 6,896 saved latent vectors
for every checkpoint and every direct OOF prediction. Maximum latent absolute
difference is 9.1553e-5; maximum direct-prediction difference is 1.2398e-5.
The CPU/GPU float32 latent comparison uses absolute tolerance 2e-4 and relative
tolerance 2e-5. The direct comparison uses 2e-5 for both tolerances. No model was
retrained or modified by this audit. LGBM fitting itself is not numerically replayed.

`8df31bae/outputs/official-metrics.csv` and `verification.json` preserve results.
Initial audit `b788441a` failed because NumPy's finite check received a Polars
Series rather than an explicit NumPy array. Audit `2095825a` reached numerical
reconstruction and tripped a stricter near-zero latent tolerance on fold 3.
Both attempts remain recorded; the successful run includes actual reconstruction
errors for all checkpoints.
