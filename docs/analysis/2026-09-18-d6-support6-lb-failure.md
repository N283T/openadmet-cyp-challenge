# CYP2D6 support6 leaderboard failure diagnosis

## Conclusion

The support6 submission failed because a very weak OOF swap advantage was
treated as if it were robust, and six highly correlated TabPFN 3.5 models gave
false confidence in a systematically different blind-test boundary. The
nearest-training support hypothesis was incomplete: overall proximity did not
distinguish the swapped sets, while proximity specifically to positive versus
negative labels strongly favored the removed control calls.

Recorded run `20260918T104416170387Z-9af2247d` joins the frozen pre/post official
leaderboards, the exact control and support6 call vectors, all six member call
vectors, OOF labels, blind-test Morgan support, standardized structures and the
canonical D6 training labels. No hidden per-compound label is inferred.

## Leaderboard movement

The fixed 49-call replacement swaps 27 control positives for 27 support6
positives. D6 MCC falls `.30864→.10973`, precision `.34121→.16545`, recall
`.38349→.19136`, F1 `.35709→.17504`, and rank `20→88`.

Using the previously fixed 7.0626% hidden-prevalence estimate, the leaderboard
bootstrap means imply a loss of approximately 10.18 true-positive equivalents
per 750 predictions and an increase of 11.91 false-positive equivalents. These
are aggregate bootstrap-weighted estimates, not recovered labels for individual
test compounds.

## The OOF selection signal was only two labels

At the fixed OOF call count, support6 and the control retain only 26 common
calls and exchange 71 in each direction. The support6-only cohort contains
31/71 positives (43.66%); the control-only cohort contains 29/71 (40.85%). The
entire apparent boundary improvement is therefore two labels among 142 swapped
OOF observations. Fisher's exact two-sided p-value is `.865`.

The highest support bin is even thinner: support6-only has 3/5 positives and
control-only 2/5. Support weighting amplified a one-observation difference in
that bin. It did not establish a stable validation advantage.

## Marginal support hid label direction

The blind swapped cohorts have essentially identical nearest-training Morgan
support: `.4854` for control-only versus `.4810` for support6-only. The
difference is `.0043`, with permutation `p=.913`. A reweighting based only on
this marginal support cannot distinguish them.

Label-aware neighborhoods do distinguish them:

| Blind-test diagnostic | Control-only | Support6-only | Control minus support6 | Permutation p |
|---|---:|---:|---:|---:|
| Nearest positive similarity | .4456 | .3600 | +.0856 | .0280 |
| Nearest negative similarity | .3743 | .4224 | -.0480 | .1568 |
| Nearest-neighbor positive rate | 74.1% | 40.7% | +33.3 pp | .0281 |
| Top-5 positive rate | 48.1% | 28.1% | +20.0 pp | .0028 |
| Top-10 positive rate | 38.5% | 27.4% | +11.1 pp | .0015 |
| Similarity-weighted top-10 positive rate | 42.2% | 28.1% | +14.1 pp | .0003 |

The 22 calls retained by both models also look label-locally positive: their
nearest-neighbor positive rate is 81.8%. The control-only calls resemble this
stable-positive set much more closely than the support6-only calls do. Gross
physicochemical descriptors are similar across the swapped cohorts; logP is the
largest difference (`3.22` versus `2.73`) but its permutation p-value is `.077`.
The failure is therefore better localized to label-neighborhood ranking than to
a broad molecular-size or polarity shift.

## Six votes were effectively one family vote

The six members are the same TabPFN 3.5-fast learner over three overlapping
feature pools and two split policies. On blind test their mean pairwise call phi
is `.865` and mean Jaccard overlap `.777`; under the exchangeable-correlation
approximation, six votes contain only about `1.13` independent votes.

Their continuous blind-test ranks are even more redundant: mean pairwise
Spearman is `.970` (range `.946–.994`). For the same feature pool, changing only
random versus scaffold gives Spearman `.994`, `.993` and `.986`. The six members
therefore differ less on deployment than their names imply.

Support6-only compounds receive 5.30 of six votes on average, while control-only
compounds receive only .44. Adding the UMAP members does not repair the ranking:
the all-nine vote retains 26/27 support6-only calls and only 1/27 control-only
calls. This is a TabPFN-family top-tail error, not primarily a random-versus-UMAP
split choice.

The member call correlation is much lower in OOF (`.623`) than on blind test
(`.865`); continuous-score Spearman likewise rises from `.856` to `.970`. The
apparent ensemble diversity collapses exactly where it is used for deployment,
so model agreement is not calibrated as confidence. Support6 and the deployed
control still have global test Spearman `.752`, but only 22/49 positive calls
overlap (call phi `.410`, Jaccard `.289`). A seemingly acceptable global
correlation therefore concealed a radically different extreme tail.

## Interpretation and decision

The blind set is closer to the labeled set in marginal Morgan similarity, but
that does not make random or scaffold validation test-like in the dimension that
matters. The missing dimension is local label composition. The control's mixed
low-capacity Morgan bank preserves positive-analogue information that the
TabPFN-family consensus systematically replaces with compounds embedded in more
negative neighborhoods.

Reject support6, all-nine TabPFN voting and nearest-any-support reweighting as D6
selection mechanisms. Keep the September 15 mean-block vector as the validated
control. For future candidates:

1. treat correlated variants of one learner as one family vote;
2. require the OOF swapped-cohort yield difference to survive uncertainty, not
   merely improve by a few labels;
3. audit nearest-positive, nearest-negative and top-k local label composition on
   the blind swaps before promotion;
4. reject large call swaps when the added cohort loses positive-neighbor evidence
   relative to the removed cohort; and
5. do not interpret marginal support matching as validation-distribution
   matching.

These diagnostics can veto an unsafe candidate without using hidden labels. They
do not reveal which individual blind compounds are positive and should not be
used to reconstruct hidden outcomes.
