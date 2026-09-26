# CYP2D6 regression OOF design after the interim result

## The design mismatch

The [official halfway post](https://openadmet.ghost.io/halfway-through-the-openadmet-cyp-challenge/)
states that the 750 test molecules came from analog expansions of strong
CYP1A2, CYP2C9 and CYP3A4 training hits. D6 hits were not an anchor. The
same 750 molecules are assayed for all four endpoints, but D6's test activity
distribution is weaker than its training distribution.

The current regression OOF is not a naive random split. Its primary validation
uses parent-retained series folds; a parent-retained scaffold split supplies a
prespecified second readout. This captures known-hit analog generalization,
and D6 was already the least stable target in the matched dual-split screen:
candidate rank correlation between the two splits was `.804`, versus
`.900/.928/.964` for 1A2/2C9/3A4. The retained anchors were already chosen
from those three other isoforms, so that aspect of test acquisition is represented.
The remaining gap is the **D6 validation population**: the observed D6 labels
near those anchors are sparse and do not match the weaker test activity mix.

## Existing other-CYP analog proxy does not solve it

Recorded coverage run `regression/20260926T004020214692Z-b2f253b0` reuses the
previous 75-anchor, 750-pair pseudo expansion selected only from measured
1A2/2C9/3A4 potency. The pairs cover 608 unique training analogs. Only
**164** have a D6 direct-inhibition label. Their mean/median D6 pIC50 are
`4.950/4.939`, versus `4.784/4.726` across all 1,493 D6 training labels.
Only 16 selected analogs are below pIC50 4. Thus this available labelled proxy
is actually somewhat **more** potent than the full D6 training set; it does
not recreate the weaker D6 test distribution described by the organizer.
Earlier series-neighbor evaluation was smaller still: 68 D6 labels, only three
below pIC50 4. A leaderboard result should not be explained or a candidate
selected using that small, opposite-direction proxy as if it were a faithful
blind-test simulation.

## Proposed validation change

Keep the existing parent-retained series OOF as a broad reference and scaffold
as the preregistered sensitivity readout. Add a **D6-specific acquisition view**
that changes the validation population, with the following constraints:

1. For each outer fold, select the 1A2/2C9/3A4 hit anchors using only the
   training side of that fold. Retain those anchors in model training, then
   identify held-out analogs by a fixed, predeclared similarity rule. Never use
   D6 outcomes to choose anchors, analogs, thresholds or weights.
2. Exclude each held-out molecule's labels from every supervised feature and
   model fit, across all heads. If the fold map changes, regenerate the nested
   log2FC features and all other fold-dependent components. Existing OOF scores
   cannot simply be reassigned to a new split.
3. Report the full D6 held-out population, the acquisition subset, and a
   separately specified low-activity stress stratum. Report label counts,
   activity distribution, ST-RAE numerator and denominator, MAE, Spearman and
   parent-group uncertainty. The low-activity stratum is an outcome-conditioned
   diagnostic, not a deployable test-selection rule.
4. Before using this view for model selection, require enough measured D6
   labels and weak-activity coverage to give stable paired comparisons. The
   current 164-label pseudo cohort fails to match the expected direction of
   shift, so it is **diagnostic only**. A genuinely dense independent D6 analog
   set would be the stronger validation source.

This is a design proposal, not a new model or a promoted selection rule. No
submitted vector, OOF fold map or candidate changes in this audit.
