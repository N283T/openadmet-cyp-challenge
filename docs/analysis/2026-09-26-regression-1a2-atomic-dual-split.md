# CYP1A2 atomic foundation complement on parent and scaffold splits

The Interim regression board puts N283T at 1A2 ST-RAE rank 32, with a gap of
`.0336` to tenth place. This audit asks whether the September 18 foundation
member has complementary signal across both existing regression splits. It
uses saved keyed OOF predictions, not a model refit or hidden-test feedback.

Recorded run `regression/20260926T015505095663Z-0acd3ef2` aligns four atomic
readouts on 1,372 nonanchor 1A2 labels: all-pool-top512 LightGBM, dense-pool-
top256 LightGBM, all-pool-top512 TabPFN 3.5-fast, and Jeremy-anchor-all
LightGBM. All are available on the same parent and scaffold fold maps. Eight
fixed standalone/equal combinations are evaluated; this is **not** a
reconstruction of the submitted ten-member ensemble.

| Readout | Parent ST-RAE | Scaffold ST-RAE |
| --- | ---: | ---: |
| Jeremy anchor LightGBM | .69999 | .69526 |
| All-pool TabPFN | .69051 | .70284 |
| All-pool LightGBM | .69337 | .69863 |
| All-pool + dense-pool LightGBM | .68780 | .69211 |
| All-pool LightGBM + TabPFN | .67947 | .68775 |
| **Jeremy anchor LightGBM + TabPFN** | **.67645** | **.67966** |

The TabPFN's standalone advantage over all-pool LightGBM on parent reverses on
scaffold, as previously reported. Its equal mean with Jeremy anchor LightGBM
is better than either atomic model on **both** splits. This supports
complementarity but does not show superiority to the full submitted bank.

Grouped paired bootstrap run `regression/20260926T015603566931Z-72c886db`
conditions on these saved predictions and resamples parent-series or scaffold
groups. For the pair minus Jeremy alone, the ST-RAE delta is `−.02355` with a
95% interval `[−.03444,−.01231]` on parent and `−.01560` with interval
`[−.02577,−.00410]` on scaffold. Against TabPFN alone, the intervals are
`[−.02657,−.00267]` and `[−.03559,−.01157]`. Against the three-model
all-pool+dense-pool+TabPFN mean, both intervals span zero. These are
post-screen, conditional intervals and exclude refit and selection uncertainty.

The incremental check `regression/20260926T015720299408Z-8e429a82`
adds either Jeremy anchor LightGBM or the complete Jeremy+TabPFN pair block
as one eleventh member to the current ten-member 1A2 mean. On parent OOF,
ST-RAE is `.670833` for current, `.670173` for the single Jeremy addition,
and `.670038` for the pair block. The pair worsens fold 0 and improves the
other four folds. The `0.000795` pooled gain is too small to warrant a new
candidate, particularly without refitting and auditing the ten-member bank
on scaffold. It also reuses an already present TabPFN member inside the new
block, so the new ensemble would weight that signal twice.

Keep the submitted 1A2 vector. The small bank increment does not justify
full-fit inference or a leaderboard probe. Preserve the dual-split atomic
finding for future bank simplification or new supervised features.
