# Classification dual-split and small-data priorities

The two user-provided small-data research reports motivated a validation audit,
a compact PLS/local-residual screen and an endpoint-specific Morgan-radius
experiment. Their recommendations were treated as hypotheses rather than as
instructions. No candidate was packaged and no submission was made.

## An exactly test-matched split is infeasible for CYP2D6

Run `5729eb3c` compares the existing random, Bemis–Murcko scaffold and UMAP
five-fold assignments with blind-test nearest-neighbor support. Morgan radius 2,
2,048-bit Tanimoto uses endpoint-labelled training molecules only. No test label
is available or used.

For CYP2D6, blind test has median support 0.471 and 42.0% below 0.4. UMAP OOF is
much harder (median 0.299; 93.4% below 0.4). Random and scaffold are closer but
remain substantially harder, with medians 0.342 and 0.342 and respectively 73.9%
and 72.3% below 0.4. Crucially, the leave-one-out training ceiling itself has
median 0.351 and 69.5% below 0.4. Therefore no honest split evaluating every D6
training row can reproduce the blind-test support distribution merely by changing
fold assignment.

For CYP3A4, scaffold is the closest actual split by Wasserstein distance (0.0468),
followed by random (0.0529), compared with UMAP (0.1695). Averaged over endpoints,
scaffold is narrowly the closest grouped choice. The existing immutable scaffold
assignment is therefore adopted as a secondary ordering check while UMAP remains
the primary OOD stress test. This is a dual readout, not permission to select the
better score after seeing results.

## PLS and local residual correction do not survive the dual readout

Run `0d2682af` evaluates five prespecified readouts on both splits: Morgan logistic
regression, raw frozen classification/mixed Chemprop embeddings, and 16-component
PLS compression of each embedding. Morgan and PLS-mixed representatives also
receive an honest local residual correction. Every outer-fold correction uses
residuals produced by inner models excluding both the outer fold and each
residual row's inner fold.

- CYP2D6 selects raw mixed Chemprop on both splits: target-prior MCC 0.0611 UMAP
  and 0.0771 scaffold.
- CYP3A4 selects raw classification Chemprop on both: 0.2114 and 0.2973.
- PLS never beats its corresponding raw embedding.
- Local residual correction slightly improves several AUROC values, but its MCC
  direction is not stable. Morgan improves on UMAP and worsens on scaffold for
  both endpoints.

Grouped bootstrap run `ff88aa1c` shows that the D6 winner is selected in only
36.1% of UMAP and 38.8% of scaffold replicates. The A4 winner is more stable at
66.4% and 64.9%. PLS compression and local residual correction are therefore
deprioritized rather than expanded into a parameter search.

## Radius-isolated unhashed Morgan is a real representation gain

Completed run `f79255ae` compares ordinary Morgan bits with endpoint-specific
radii derived from the preceding attribution analysis: radius 2 for D6 and radii
0/1 for A4. Exact sparse Morgan environment IDs use a fold-fitted `DictVectorizer`;
held-out and test-only environments are ignored rather than hashed into occupied
positions. Failed precursor `5262c679` records a scikit-learn sparse-index dtype
error; the completed run changes only those matrix indices to 32-bit integers.

For D6, the fixed rank mean of ordinary Morgan and unhashed radius-2 counts moves
target-prior MCC from 0.0384 to 0.0999 on UMAP and from 0.0677 to 0.1066 on
scaffold. The direction also holds for AUROC: 0.5729 to 0.6077 and 0.6154 to
0.6493. In grouped selection bootstrap `dc3d7523`, ordinary Morgan wins only 0.46%
of UMAP and 0.24% of scaffold replicates; radius-specific variants take the rest.
This supports the representation family even though the exact best radius recipe
varies.

For A4, ordinary-plus-unhashed radii 0/1 improves AUROC on both splits and is best
on scaffold by target-prior MCC. Unhashed-only is best on UMAP. This is useful
feature diversity, but it is not compared here against the stronger complete A4
ensemble.

## Comparison with the live D6 ensemble

Post-screen run `29814833` compares fixed late-rank mixtures with the live D6
rank. The new radius rank alone scores projected MCC 0.09988 versus 0.09485 for
live. All three mixtures are worse than live. The standalone delta is only
+0.00503, with a conditional 10,000-draw UMAP-group interval of
[-0.04135, 0.05579]. Thus the new feature family is validated against ordinary
Morgan across both splits, but replacement of the current heterogeneous ensemble
is not established.

## Fixed radius–Chemprop block fusion

Run `1bdad261` combines the validated radius block with the best raw Chemprop
block using fixed rank weights. D6 radius remains best on UMAP/scaffold
(0.09988/0.10657); every Chemprop mixture is worse. For A4, the 1:1
radius/classification-Chemprop mean improves both splits to 0.25522/0.35616,
compared with 0.22396/0.29571 for radius and 0.21142/0.29726 for Chemprop.
Grouped-bootstrap run `f88cc4f6` supports the A4 fusion family, although UMAP
does not uniquely identify the exact weight.

Run `63db076c` compares that A4 block with the complete live Caruana rank. Even
the lightest 3:1 Caruana/block mixture lowers target-prior MCC by 0.00772 and
swaps 14 calls each way. The block is therefore an internal baseline, not a
live-model complement.

## Primary FDA/Instem CYP3A4 TDI data

The primary Faramarzi et al. supplement was recovered from Figshare after
rejecting two misleading routes. Direct PMC run `b73557c5` received an HTML
challenge page. Run `a0255df8` downloaded an unrelated pharmacovigilance archive
from a mismatched BioStudies `DataSheet1` path and is unusable. Correct fetch run
`594c2424` enforces the article XML's 3,331,650-byte size and MD5
`995f8effde82473e0c5625df1ebe9420` before extraction. Nested-launch failure
`94e6b294` is retained.

Audit `47a95438` finds 621 unique binary CYP3A4 TDI connectivities, 49.28%
positive. Forty exactly overlap challenge structures, all on the train side and
none in blind test. Of 17 overlaps with an A4 challenge label, only 11 agree
(64.7%). This supports assay/domain shift and rules out naïve label pooling. The
one-record difference from the paper's stated 620 compounds is preserved rather
than silently edited. Failed audits `c96d0012` and `65a3776a` record a missing
derived-connectivity assumption and a non-UTF-8 metadata field.

Run `5463a4dc` excludes all 40 overlaps and trains LR, ExtraTrees and a
Tanimoto-kernel SVM on the remaining 581 source molecules as an independent
task. The best UMAP mixture reaches 0.27860 versus 0.25522 for the challenge
block, but the same LR mixture falls to 0.32210 versus 0.35616 on scaffold.
Selection bootstrap `4ec036d6` gives the UMAP winner only 27.5% frequency;
scaffold retains the challenge-only block in 70.1%.

Similarity-gate run `d3b629c5` shows low source support (median maximum Morgan
Tanimoto 0.250 overall and 0.260 in blind test). Restricting FDA influence to
higher-support molecules gives tiny, directionally consistent AUROC/AP gains,
but changes no fixed-rate calls and leaves MCC unchanged. Retain the source and
predictions for later representation learning; do not alter live A4.

## Alert-centered local environments

Run `b95e8380` isolates collision-free Morgan environments whose centers lie on
the 12 explicit TDI-alert matches or their one-bond neighborhoods. D6 improves by
only 0.00510 on scaffold and worsens on UMAP; A4 mixtures worsen on both splits.
The direct alerts therefore do not become useful merely by replacing counts with
local exact environments. A true metabolism-site model remains a distinct,
untested hypothesis; this approximation is rejected.

## Decision

Use UMAP plus scaffold for subsequent low-cost comparisons and report both. Stop
PLS, hard/local support correction and explicit-alert local-environment work.
Keep the A4 radius–Chemprop block as an internal baseline and the FDA source as
domain-shifted auxiliary data, but neither complements live A4. Promote
endpoint-specific unhashed Morgan only as a future D6 stack block. No candidate
was packaged and no submission was made. A true atom-level metabolism-site model
remains distinct from the rejected alert approximation.
