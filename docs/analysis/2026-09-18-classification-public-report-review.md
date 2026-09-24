# Classification public-report review

## Decision

Freeze the scored CYP3A4 vector and resume classification work on CYP2D6. The
three public reports named for this review all beat the current CYP2D6 score,
but most of their broad ideas are already represented in this project's saved
experiments. The one bounded, materially untested recipe is vaishnavi53's
shallow Random Forest / Extra Trees / LightGBM probability mean. Test that
architecture as a new keyed OOF peer, with and without the existing approximate
TDI-alert panel, before doing any further representation or threshold search.

Do not reopen global threshold tuning, direct alert concatenation by itself,
TabICL on Jeremy600, or continuous-shift regression. Those routes are already
covered by project evidence, and the public reports do not overturn it.

This review performs no model fitting, candidate creation or submission.

## Sources and current board positions

The source snapshot is the dashboard generated on 2026-09-18. Scores are the
full values already retained in the corresponding local leaderboard captures;
the table displays rounded values for readability.

| Participant | CYP2D6 MCC/rank | CYP3A4 MCC/rank | Macro MCC/rank | Public source |
|---|---:|---:|---:|---|
| N283T | .3086 / 19 | .5299 / 1 | .4193 / 3 | current project |
| vaishnavi53 | .3412 / 9 | .4679 / 6 | .4045 / 7 | Google methodology report |
| 450nm | .3408 / 10 | .4775 / 5 | .4091 / 5 | campaign console |
| jeremy | .3401 / 12 | .4648 / 8 | .4025 / 8 | public reproduction repository |

All three reports are relevant because they beat N283T on CYP2D6, not because
they beat N283T overall. The D6 gaps are only about .0315--.0325 MCC and each
entry's marginal bootstrap standard deviation is about .085--.091. These are
not paired uncertainties, so the board does not establish that one named recipe
is intrinsically better. It does establish that D6 remains the only useful
endpoint to move while A4 is rank 1.

Archived sources:

- `references/community/cyp/450nm-campaign-console/20260918T051404Z/`
- `references/community/cyp/vaishnavi53-methodology/20260918T051404Z/`
- `references/community/cyp/jeremy-reproduction-kit/20260918T051404Z/`

Jeremy's local clone was already present at commit
`d3a5549e0c98167435cbf8f5115c68a9a184ed7a`; the archive records hashes for the
README and relevant TDI scripts without duplicating the unlicensed source code.

## vaishnavi53

The reported D6 pipeline uses ECFP4 (radius 2, 2,048 bits), 13 binary
mechanism-based structural alerts, five Butina folds, and an equal probability
mean of:

- Random Forest: 600 trees, maximum depth 8;
- Extra Trees: 600 trees, maximum depth 10;
- LightGBM: 300 estimators, maximum depth 4, 15 leaves, learning rate .05.

The final decision threshold is selected from OOF predictions for MCC. The
report does not disclose the 13 exact SMARTS, fold assignment, random seeds,
class-weight settings for D6, threshold grid, or whether the descriptors
described earlier are also present in the D6 branch. Exact reproduction is
therefore impossible from the report alone.

The project already tested a transparent 12-category alert approximation.
Several individual motifs were enriched, especially methylenedioxyphenyl on D6,
but concatenation worsened honest inner-threshold macro MCC. Extra Trees was
also screened as a standalone fingerprint family and had useful AUROC but weak
target-prior MCC. What has not been tested is the report's deliberately shallow
three-family probability mean. Its potential value is model-family smoothing,
not evidence that the alert block itself works.

## 450nm

The most important public result is operating-point diagnosis rather than a new
representation. The report recovers the scored CYP2D6 prevalence as 7.063% from
published accuracy/precision/recall/MCC and observes that the then-current top
ten entries called 2.4--8.3% positive, median 6.54%. Moving their D6 call rate
from about 21.6% to 7% raised MCC from .1862 to .3408 without changing the model.

This project independently inferred the same 7.06% hidden prevalence and its
current D6 vector already calls 49/750 = 6.53% positive. There is no unclaimed
threshold gain here. In particular, the report is evidence against another OOF
threshold search: its own nested threshold optimization reduced D6 MCC, and its
local ROC recommendation pointed in the wrong direction on the board.

The report also finds that the continuous TDI shift is noise dominated. On
CYP2D6, its reported scaffold-CV Spearman is .076 even though the two component
pIC50 levels are more predictable; shift regression, uncertainty weighting and
separate-level subtraction all lose to direct binary classification. Treat this
as external negative evidence against opening a continuous-shift modelling
branch for D6 classification.

Direct TDI-alert concatenation is likewise negative in that campaign. The
author reports motif enrichment but no stable MCC gain, consistent with this
project's controlled alert screen.

## jeremy

The public kit's TDI branch is LightGBM on ECFP4 plus RDKit descriptors and
TabICL8 on a frozen 600-dimensional `chemprop_medium` embedding reduced to 200
components. It uses common Butina folds, Caruana bagged ensemble selection and
an OOF-selected MCC threshold.

This is already deeply represented here:

- the Jeremy600 embedding is archived and has matched TabICL8/TabPFN8 reads;
- Jeremy600 TabICL is the strongest tested TabICL representation by D6 AUROC;
- it is already present in the broader OOF prediction bank and historical A4
  ensemble work;
- multiple Jeremy-derived and project-owned replacement encoders have been
  compared, distilled and adapted.

There is no missing broad Jeremy experiment. The repository's tracked
`results/submission_tdi.csv` contains 29 D6 and 216 A4 positive calls, whereas
the public leaderboard row's precision and recall imply that the submitted D6
vector was not that exact tracked file. The file is therefore not labeled as a
leaderboard-proven prediction vector and should not be blended on that premise.

## Next experiment

Implement one recorded D6-only comparison against the current classification
OOF bank:

1. use the canonical keyed classification labels and fixed UMAP folds, with the
   established scaffold split as a secondary ordering check because local/LB
   divergence is already material;
2. fit the three reported shallow families on ECFP4, first without alerts and
   then with the existing documented approximate alert panel;
3. save each atomic OOF/test probability and the fixed equal mean; do not tune
   ensemble weights;
4. compare AUROC, AP and fixed-rate MCC at 49/750-equivalent prevalence, not an
   OOF argmax threshold;
5. add the equal-mean peer to the current D6 prediction bank only if it improves
   both primary and secondary ordering or provides a clearly complementary
   fixed-rate boundary under grouped uncertainty.

This isolates the only new hypothesis in the reports: a shallow heterogeneous
tree mean may transfer better than either a single alert-augmented learner or a
larger selected ensemble. It also leaves the rank-1 A4 vector untouched.

## Experimental follow-up

Runs `71dcd147`, `02baa35e`, `8dae5f36` and `31129d64` execute the bounded
follow-up. CYP3A4 is untouched throughout.

The vaishnavi53-style equal probability mean is rejected. On the primary UMAP
split its best equal-three variant reaches target-prior MCC `.07526`, below the
current `.09988` control used in that run. Adding the 13 approximate alert
features does not rescue it. The best alert-augmented atomic model is the
balanced Random Forest: its UMAP AUROC `.58683` exceeds the control `.57926`,
but its fixed-rate MCC is only `.05188`. On scaffold it instead reaches
MCC/AUROC `.13265/.62381`, a large split reversal rather than dual-split
confirmation.

Run `02baa35e` then transfers the strongest recent regression feature pools into
fold-safe D6 LightGBM classifiers. `anchor_all_pub030_raw::balanced` gives the
best primary ranking, AUROC `.59630` versus `.57926` for the control, but its
fixed-rate MCC is `.08008`. `owned_scratch_anchor::balanced` is strong only on
scaffold, MCC/AUROC `.12732/.63915`. These are complementary ranking signals,
not standalone replacements.

Run `8dae5f36` reconstructs the currently scored mean-block D6 OOF recipe and
tests three prespecified peers at candidate rank weights 10%, 20%, 25%, one
third and 50%. Only the 10% `anchor_all_pub030_raw::balanced` addition improves
the point estimates together:

| Model | AUROC | AP | target-prior MCC | blind-test swaps each way |
|---|---:|---:|---:|---:|
| Deployed mean-block control | .58094 | .29506 | .11011 | 0 |
| + 10% anchor/pub030 regression rank | **.58627** | **.29512** | **.11531** | 6 |

The grouped MCC delta is only `.00527 [-.01762,.02876]`, with 66.3%
improvement probability. Fold-excluded weight selection in `31129d64` chooses
the 10% blend for held folds 0, 3 and 4 and the control for folds 1 and 2. Its
honest combined AUROC/AP improve to `.58277/.30023`, but target-prior MCC falls
to `.10496`. Thus the regression rank has repeatable broad-ranking value but
does not stably improve the 49-call boundary.

No classification candidate or submission is created. Preserve the deployed
D6 and A4 vectors. The next D6 experiment should introduce a genuinely new
source of matched TDI-shift supervision or a new model family, rather than
further manual weighting of these public-report and regression-feature peers.
Continuous direct/TDI level prediction, metabolite activity, auxiliary
regression fine-tuning and call-rate correction are already covered by earlier
negative or non-incremental experiments.

## Newly released TabPFN 3.5 follow-up

The preceding classification bank used TabPFN 2.6. After the user identified
the newly released TabPFN 3.5 as a distinct untested model family, runs
`e0583fe4` and `387d4e0e` evaluate the local 3.5-fast and full classifier
weights with eight estimators. Five regression-derived feature pools are
screened with 3.5-fast; the best three are confirmed with full 3.5. LightGBM
top-512 feature selection and missing-value imputation are repeated inside each
outer training fold.

TabPFN 3.5 produces the strongest new D6 ranking in this follow-up:

| Backend / features | UMAP AUROC | UMAP AP | UMAP target MCC | scaffold target MCC |
|---|---:|---:|---:|---:|
| Current fixed-ratio control | .57926 | .29124 | **.09988** | -- |
| 3.5-fast, Jeremy + anchor | .60386 | .31538 | .09485 | .09654 |
| 3.5-fast, anchor + pub030 | **.61262** | .31332 | .08987 | .12732 |
| 3.5-fast, owned + scratch + anchor | .61179 | **.31615** | .07049 | **.15465** |
| Full 3.5, Jeremy + anchor | .60448 | .30866 | .09485 | .12732 |
| Full 3.5, anchor + pub030 | **.61365** | **.31571** | .08495 | .12204 |

Thus the new model family improves broad ranking on both structural views, but
no standalone readout beats the control at the UMAP 49-call boundary. Fixed
blend audits `c1058115` and `d40335d6` show that the most useful direct addition
is 20% 3.5-fast Jeremy+anchor: it preserves deployed mean-block MCC `.11011`
while improving AUROC/AP `.58094/.29506→.59021/.30946` and swaps five blind
calls each way. The full model does not improve this boundary.

Run `bc666ebe` averages the matched fast and full Jeremy+anchor ranks. A 25%
addition again preserves MCC `.11011` and improves AUROC/AP to
`.59238/.31018`, with eight blind swaps each way. The two 3.5 backends are
highly redundant (OOF Spearman `.97142`), and grouped MCC delta remains centered
near zero (`-.00073 [-.02760,.02562]`). These percentage sweeps are diagnostic
sensitivity checks, not research candidates: promoting them would violate the
existing OOF-bank rule against inventing a manual blend location for each new
expert. Add the TabPFN 3.5 scores to the keyed bank for a prespecified bank-level
comparison instead. No candidate or submission is created and A4 remains
unchanged.

## CYP3A4 TabPFN 3.5 symmetry check

At the user's request, run `c6b07258` applies the same 3.5-fast protocol to
CYP3A4 despite its current leaderboard rank 1. The check is informative rather
than a reason to relax the endpoint freeze. The all-feature top-512 pool is
strong on both splits:

| Model | UMAP AUROC | UMAP AP | UMAP target MCC | scaffold target MCC |
|---|---:|---:|---:|---:|
| Development fixed-ratio control | .80713 | .51273 | .45467 | -- |
| 3.5-fast all-pool | **.81689** | **.55173** | **.46999** | **.50516** |
| 3.5-fast dense-pool | .81477 | .55125 | .45007 | .50058 |

The relevant scored control is the later equal-seven A4 model, not the earlier
fixed-ratio development control. Audit `6628e253` therefore compares directly
with its exact OOF and test score. A one-third 3.5-fast all-pool rank addition
raises AUROC/AP/MCC `.80571/.52160/.46540→.81617/.54550/.47458`, swapping 13
of 249 blind calls. The grouped MCC delta is `.00898 [-.00349,.02103]`, with
92.2% improvement probability.

Full 3.5 confirmation `5c3e72b5` gives all-pool UMAP
AUROC/AP/MCC `.81893/.55017/.46080` and scaffold MCC `.50516`. Correctly
labeled equal-seven blend audit `ce616a0d` reproduces the one-third improvement
at `.81657/.54676/.47305`, with 14 blind swaps and 87.1% grouped improvement
probability (`.00753 [-.00567,.02017]`). Superseded audit `bcab335e` has the
same values but inherited a fast-model output label from the generic script.

This is a real, backend-replicated A4 ranking signal, but neither conditional
interval excludes zero and the current A4 leaderboard gain already greatly
exceeded its CV projection. Preserve the rank-1 A4 submission vector. Retain
the one-third result only as evidence that TabPFN 3.5 contains complementary
signal. The manual blend is not a future candidate; any reconsideration must be
through a prespecified common-bank procedure. No candidate, submission or
upload is created.

## No-manual-weight bank follow-up

Run `f7c4b04c` appends the five 3.5-fast feature-pool predictions to the common
D6 keyed OOF bank and compares the legacy and extended banks under identical
equal-weight procedures. The all-expert mean improves AUROC
`.58323→.59169` but reduces target-prior MCC `.09988→.09485`. Nested
outer-fold-excluded equal-weight subset selection also worsens MCC
`.09485→.08008`; four of five outer selections contain no new 3.5 expert.
Against the scored mean-block control `.11011`, the extended nested subset has
grouped delta `-.02988 [-.06848,.00518]`. This rejects the bank addition without
introducing a manual blend percentage.

Standalone run `19e2a462` increases only the TabPFN 3.5-fast internal estimator
count from 8 to 32. All three UMAP top-512 candidates lose fixed-rate MCC; the
best is `.08495`, below both the corresponding 8-estimator models and the
control. Minor scaffold gains reinforce split instability rather than a robust
tail improvement. TabPFN 3.5 remains useful ranking evidence, but neither bank
addition nor estimator stabilization produces a D6 candidate.

## Test-like validation follow-up

The Morgan representation family was already broadly audited; the material
unresolved difference is validation support. D6 median nearest-training Morgan
similarity is `.299` for UMAP, `.342` for random folds and `.471` for blind test.
Random-fold refits reverse the UMAP conclusion: TabPFN 3.5-fast anchor+pub030
scores target-prior MCC `.12590`, Jeremy+anchor `.12058`, while the public
shallow-tree equal mean reaches `.11011` (`8b586046`, `3f5b3fc3`).

Support post-stratification `892612f4` ranks Jeremy+anchor first at weighted MCC
`.14874`, but the estimate relies on only 35 random-OOF rows above similarity
`.6` to represent 20.1% of blind test. Equal split-bag audit `fb827f03` finds
scaffold all-pool strongest by ordinary MCC `.13675` and three-split Jeremy
stronger under support weighting at `.16101`. Complete-grid top-call consensus
`8829d528` improves the scored mean-block from ordinary/support-weighted MCC
`.11011/.10398` to `.12058/.12090`, with 27 of 49 blind calls changed.

The evidence supports validation mismatch—not an untried Morgan variant—as a
plausible explanation for the public leaderboard gap. These are high-change
prospective LB probes rather than a selected candidate; no package or upload is
created.

## New random split search

Support-only searches over 500 D6-specific seeds find ordinary KFold seed 215
closest to blind test at Wasserstein distance `.09681`, but its fold positive
rates vary from 17.4% to 24.0%. StratifiedKFold seed 498 is nearly as close at
`.09760`, raises median OOF support to `.34483`, and fixes every fold near the
overall 21.6% positive rate (`0e8db36d`, `2fe346d6`). Blind-test median support
is still `.47114`, so split search cannot reproduce the deployment geometry.

Matched TabPFN refits across the five support-selected stratified seeds show
large residual variance: anchor+pub030 MCC ranges `.08519–.14688` and
Jeremy+anchor `.09002–.12538` (`277c5b65`, `03186dea`). Exact global-boundary
audit `83ba948f` averages all five label-blind selected seeds. The anchor bag
only moves MCC `.10832→.11347`; the Jeremy bag ties `.10832`, and each changes
26 of 49 blind calls. Seed 498's strong standalone point estimate does not
survive seed bagging. Use this split family for secondary stability checks, not
best-seed selection; no candidate or upload is created.

## Submission disposition

The complete six-model random-plus-scaffold top-call vote was selected as the
cleanest prospective D6 leaderboard probe because it improves both ordinary and
support-weighted validation without a hand-tuned score percentage. Packaging
run `d535a67b` produced 49 D6 calls and changed 27 calls in each direction from
the scored vector. The companion A4 prediction was copied bitwise unchanged.

Candidate `20260918T094447504193Z-d15100fb`, prepared submission
`20260918T094451131332Z-f23dd61b`, and live attempt
`20260918T094516432346Z-291aa76a` record the submission chain. The API accepted
the file with status `sent`; the result remains pending leaderboard processing.

## Leaderboard result

The probe fails decisively. D6 MCC drops from `.30864` to `.10973` and rank
from 20 to 88; precision, recall and F1 all roughly halve. Macro MCC drops from
`.41927` to `.31982` and rank from 5 to 38. A4 remains exactly `.52991`, rank 1,
which confirms that the package isolation worked and attributes the degradation
to the D6 replacement.

Reject the support6 consensus. Its random/scaffold ordinary and support-weighted
validation gains are not predictive of this blind operating tail. The September
15 mean-block D6 vector remains the leaderboard-validated control; no further
submission should be selected from this split-consensus family without new,
independent evidence.
