# Regression dual-split feature audit

## Why add a second split now

Routine regression development previously used only the parent-retained series
split. That was an explicit early decision: reconsider split sensitivity only
after leaderboard and CV diverged materially. The latest submission supplied
that evidence. Its new CYP2D6 bank improved OOF but worsened leaderboard ST-RAE
from `.572470` to `.582788`; CYP1A2 also worsened while CYP2C9 and CYP3A4
improved.

The classification/TDI dual-split work had already established the intended
role of the second readout. Scaffold was the grouped split closest to blinded
test structural support on average, while UMAP was a harder OOD stress test.
The two scores were reported together; the better score was never selected
post hoc. Regression now adopts the same principle, adapted to its current
membership design:

- parent-retained series remains the primary development split;
- parent-retained scaffold is a prespecified secondary ordering check;
- improvements and candidate ranks are reported target by target on both;
- the scores are not averaged, and a candidate is not chosen because it wins
  only the more favorable split.

## Secondary membership

Run `11c1a26e` joins the immutable scaffold assignment from `edba6fc0` to the
same parent metadata used by primary run `a3923255`. All 70 anchor compounds
remain in training in every outer fold. The assignment covers all 6,146 train
compounds, contains five folds and has zero connectivity groups crossing folds.

## Matched LightGBM screen

Run `c89b25fa` refits 33 feature candidates from primary screen `5c46e72d` on
the secondary scaffold folds. LightGBM parameters, target labels, feature
blocks and outer-training-only feature selection are unchanged. Four candidates
using nested predicted log2FC are excluded because those features were trained
for the primary folds and cannot be reused honestly on scaffold. Run `ef124be7`
performs the matched ordering audit.

Lower ST-RAE is better.

| Target | Rank correlation | Primary winner | Series | Scaffold rank / score | Scaffold winner | Scaffold | Series rank / score |
|---|---:|---|---:|---:|---|---:|---:|
| 1A2 | .900 | all_pool_top512 | .693368 | 3 / .698627 | jeremy_anchor_all | .695263 | 5 / .699994 |
| 2C9 | .928 | owned_scratch_anchor | .489581 | 5 / .502317 | dense_pool_top512 | .491042 | 3 / .496386 |
| 2D6 | .804 | anchor_all_pub030_raw | .832272 | 13 / .837916 | jeremy_anchor_all | .822203 | 3 / .832777 |
| 3A4 | .964 | dense_pool_top256 | .422655 | 4 / .414048 | all_pool_top256 | .411148 | 3 / .427181 |

The broad conclusions are more stable than the exact winners. The four selected
dense pools occupy ranks 1-4 on both splits for 3A4. For 1A2, the primary winner
remains third on scaffold, and `dense_pool_top256` is second on both. For 2C9,
the primary winner falls to fifth but remains in the leading group; the two
dense pools are more balanced across splits. CYP2D6 is the least stable target:
its primary non-log2FC winner falls from first to thirteenth, while
`jeremy_anchor_all` ranks third and first. This agrees with the observed D6
OOF-to-leaderboard failure and argues against another single-split D6 replacement.

Relative to the fixed Morgan-plus-Mordred control, 29/33 candidates improve on
both splits for 1A2 and 30/33 do so for each other target. Thus the new checkpoint
features are not merely an artifact of the primary split, but fine-grained
winner selection is split-sensitive.

## Decision

Use the two regression readouts prospectively. For low-cost screens, require
both target-specific score tables and prefer candidates that remain in the
leading group on both splits. Do not require the exact same winner, average the
two scores, or turn scaffold into a post-hoc fallback. Before changing the
production banks, regenerate fold-dependent nested log2FC on scaffold and then
evaluate the reduced target-specific bank candidates on both splits. Prioritize
CYP2D6 because it has the lowest rank agreement and a confirmed leaderboard
transfer failure.

