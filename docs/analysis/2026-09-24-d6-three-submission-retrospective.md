# September 15 D6 versus support6 and AP-mass90

## Scope and evidence

Recorded run `classification/20260923T190811462293Z-4f4a92f2` joins the exact
accepted September 15, September 18 support6, and September 23 AP-mass90
submission vectors with their official D6 leaderboard snapshots. It also joins
the saved honest OOF decisions and the previously computed blind-test Morgan
label-neighborhood diagnostics. A second recorded run,
`classification/20260923T191425347911Z-9a9f9ab2`, tests a separately trained,
fold-excluded primary-screen response predictor as an auxiliary explanation.
The two failed attempts `9758b270` and `2e88b08a` produced no scientific
result; they document correction of a split-ID join and an integer-label type.
Completed run `9ec3ebc2` was superseded after its anchor exception was found.

The three packages each call 49 of 750 compounds positive for CYP2D6 and have
identical CYP3A4 calls. Thus their D6 leaderboard differences come from the
identity of the called compounds at this operating point. This does not prove
that 49 is the best call count.

| Submission | D6 MCC | D6 precision | D6 recall | Difference from September 15 |
|---|---:|---:|---:|---:|
| September 15 control | .30864 | .34121 | .38349 | reference |
| September 18 support6 | .10973 | .16545 | .19136 | −.19891 |
| September 23 AP-mass90 | .23483 | .27996 | .30522 | −.07381 |

Leaderboard metrics are bootstrap means on an unknown confidently labelled
subset. They do not identify individual hidden labels. Under the previously
estimated 7.06% hidden prevalence, the recall differences correspond roughly to
10.2 and 4.1 fewer true-positive equivalents per 750 for support6 and AP-mass90,
respectively; these are *not* recovered true-positive counts.

## Exactly which calls changed

Support6 removed 27 September 15 calls and added 27. AP-mass90 removed 15 and
added 15. All 15 AP-mass90 removals are among support6's 27 removals. AP-mass90
kept the other 12 control calls and avoided 21 of support6's new calls. This
explains why AP-mass90 is less damaged than support6, but does not alone explain
why September 15 still wins.

On the honest labelled-only OOF decisions, the support6 exchange is 71 calls
each way: the additions contain 31 positives and the removals 29 (Fisher
two-sided `p=.865`). AP-mass90 exchanges 48 each way: 18 positives added and 17
removed (`p=1.0`). Both apparent improvements are only one or two labels among
many swapped observations. This comparison is a retrospective fixed-call
diagnostic; the OOF recipes differ from their full-fit deployment in some
folds, and the labelled-only population is only 1,497 of 6,146 training
molecules. It cannot reproduce the blind top-49 selection.

## Why support6 failed; why that diagnosis does not cover AP-mass90

On blind-test compounds support6 removed by itself, the similarity-weighted
positive rate among the ten nearest labelled Morgan neighbors is `.422`. Its
added compounds average `.281` (two-sided permutation `p=.0006`). Thus the
September 15 calls retain substantially stronger local positive-label evidence.
Six support6 members were also near-duplicates on the blind set (mean pairwise
rank correlation `.970`); their agreement overstated independent support.
These observations give a plausible, specific account of the support6 loss,
though hidden labels prevent causal proof.

For AP-mass90 the same Morgan statistic is `.37509` for the 15 removed control
calls and `.37513` for the 15 additions (`p=.9992`). Nearest-positive similarity
actually favors the AP additions (`.514` versus `.380`, `p=.0322`). Therefore
the neighborhood rule that identifies support6's mistake **would have favored
or accepted AP-mass90**, and cannot explain its leaderboard loss. Retrospectively
turning that rule into a universal candidate gate would repeat the evaluation
problem.

## Auxiliary primary-screen proxy check

The saved CYP2D6 primary log2FC predictor was fit without each nonanchor
training compound's own UMAP fold. Its assignment differs from the
classification OOF assignment, so the audit joins by compound ID and uses the
primary model's own excluded fold. The 36 labelled anchor compounds remain
eligible for primary-model training and are excluded from this proxy's OOF
accuracy and contrast calculations. For blind compounds it averages all five
outer predictions. Lower log2FC predicts stronger primary inhibition.

It is **not** an accuracy proxy for the D6 TDI class. Against the supplied D6
OOF class labels for the 1,461 nonanchor molecules, `−predicted log2FC` has
AUROC `.443`; class positives average `−.890`, negatives `−1.052`. Strong
primary inhibition is not the same event as
time-dependent inhibition among the labelled molecules. The primary predictor
also appears in feature pipelines used by some classification members, so it is
not a wholly independent model-family vote.

The September 15-only calls have lower predicted primary log2FC than the
support6-only and AP-mass90-only calls both in OOF and blind test. For the AP
blind swap the difference is `−.112` with permutation `p=.639`; on nonanchor
OOF it is `−.347`, `p=.0178`. Given the reversed OOF class association and weak blind
contrast, this proxy does not explain why September 15 beats AP-mass90. The
negative result is retained to prevent using primary inhibition as a
convenient but unsupported selection signal.

## Decision and next evaluation

The retrospective evidence supports a support6-specific failure mechanism and
shows exactly what AP-mass90 fixed relative to support6. It does **not** isolate
the molecular or assay reason for AP-mass90's remaining loss to September 15.
Its OOF gain was tiny and uncertain, its blind additions passed the known
label-neighborhood check, and its stable outer-fold votes established
repeatability rather than correctness. A single aggregate leaderboard result
cannot assign the AP errors to specific compounds or assay regimes.

The next candidate evaluation must first refit the September 15 and AP-mass90
recipes on the **same alternative outer split**, preferably one chosen before
inspecting either model's result and closer to the blind acquisition geometry.
Existing UMAP OOF predictions cannot be relabelled as scaffold or random OOF:
their training sets would include many of the newly held-out molecules. This
experiment directly tests whether September 15's advantage appears when the
models are evaluated away from the UMAP split that preferred AP-mass90. The
available scaffold scores cover some AP-mass90 member families, but not the
complete AP recipe or the deployed September 15 hierarchy.

Then generate outer-fold scores for the complete 6,146 training-molecule pool,
keeping null labels for unscored compounds; select calls within that full
held-out pool and calculate MCC only on scored rows. Compare candidate and
September 15 control at jointly selected call rates, with paired swap yields
and uncertainty on multiple prespecified split families. Inspect
whether the added blind calls have supporting label neighborhoods and whether
model-family agreement collapses on blind data. A candidate with a plausible
gain and wide uncertainty can be a deliberate leaderboard probe; these
diagnostics should quantify the risk and predicted gain rather than form an
automatic veto. Neither the exact alternative-split recipe comparison nor the
complete-pool backtest has yet been run.
