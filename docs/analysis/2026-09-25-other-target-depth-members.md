# Other regression targets: intermediate-depth members in the submitted ensembles

## Question and design

After the CYP2D6 whole-member screen, test the same intermediate CheMeleon and
Jeremy FT600 ideas on CYP1A2, CYP2C9 and CYP3A4. The comparison uses the actual
submitted endpoint recipe and matched parent-series five-fold OOF. No fixed
fractional blend is searched. A CheMeleon depth block is the equal mean of its
depth-specific LightGBM predictions **inside one candidate member**. The Jeremy
block is formed in the same way. Outer additions and replacements operate on
whole members.

The three submitted recipes require different treatment:

| Target | Submitted recipe | OOF control |
| --- | --- | ---: |
| 1A2 | Equal mean of ten atomic members | .670833 |
| 2C9 | Five-member equal subset selected from eleven candidates on each outer training fold | .454711 |
| 3A4 | Caruana selection from seventeen candidates, with steps selected by inner CV | .394567 |

The exact saved OOF for each control was reconstructed from its atomic inputs
and fold-specific fits to less than `1e-9` maximum prediction difference.
For 1A2 and 2C9, fixed additions enter as one new equal member. Fixed
replacements preserve the incumbent's weight. For 3A4, fixed replacements
preserve its saved per-fold Caruana weight. Runs `04d782a3` and `b931aa2a`
also evaluate 2C9/3A4 by refitting their original selection rules after a
candidate is added to the bank or substituted for its corresponding member.

CheMeleon depth-1 TabPFN and TabICL use fold-training PCA256 and nested
log2FC4 (`c352d23e`, `a53c9eef`). The matched CheMeleon depth 1/3/6
LightGBM OOF is `dcf9b1b8`. Jeremy FT600 depth 1/2/4 LightGBM readouts use
the actual fold-specific submitted encoders and the submitted nested log2FC4
recipe (`22c133cf`, `b72b5dc0`, `d9bc0461`). The final-depth control uses the
exact saved feature files and reproduces the original FT member exactly.

## Results

Lower ST-RAE is better. Jeremy FT intermediate depths are weaker than its
final depth on all three targets when used alone:

| Target | Depth 1 | Depth 2 | Final depth 4 |
| --- | ---: | ---: | ---: |
| 1A2 | .773005 | .751995 | .720947 |
| 2C9 | .580622 | .568074 | .494596 |
| 3A4 | .492291 | .479532 | .430372 |

The target-specific outer ensemble comparisons are:

| Target | Comparison | OOF ST-RAE | Change from control |
| --- | --- | ---: | ---: |
| 1A2 | Current ten | .670833 | — |
| 1A2 | Replace CheMeleon PCA TabICL with depth 1+3+6 LightGBM block | .670757 | −.000076 |
| 1A2 | Add the depth 1+3+6 block as eleventh member | .670760 | −.000073 |
| 1A2 | Replace Jeremy FT with depth 1+2+4 block | .671542 | +.000709 |
| 2C9 | Current selected five | .454711 | — |
| 2C9 | Best fixed addition, Jeremy depth 1+4 block | .456052 | +.001341 |
| 2C9 | Add any of three depth candidates and reselect the equal subset | .454711 | zero; all rejected |
| 2C9 | Replace current FT with depth 1+2+4, then reselect | .461257 | +.006545 |
| 3A4 | Current Caruana bank | .394567 | — |
| 3A4 | Best fixed replacement, Jeremy depth 1+4 block | .395173 | +.000606 |
| 3A4 | Add Jeremy depth 1+2+4, then reselect Caruana | .394822 | +.000255 |
| 3A4 | Add CheMeleon depth 1 TabICL or 1+3+6 LightGBM, then reselect | .394567 | zero; both rejected |

The 2C9 equal-subset selection never picks any of the three added depth
candidates across five outer folds. The 3A4 Caruana fit picks the added Jeremy
block in one fold, but OOF ST-RAE still rises. It never picks either added
CheMeleon candidate. Replacing the corresponding incumbent and rerunning
selection worsens both targets. The full scenario tables and fold fits are in
`artifacts/tracked/other-regression-depth-members-20260925/`.

For 1A2, parent-series-group bootstrap `ec41c034` (2,000 replicates,
1,004 groups) gives the depth 1+3+6 **replacement** a paired delta interval
`[-.004247, +.003872]` and improvement in 51.35% of resamples. The
**addition** interval is `[-.003066, +.002977]`, improvement in 50.85%.
These tiny observed gains are consistent with noise.

## Decision

Keep all three submitted regression endpoint recipes. Intermediate Jeremy
depths do not displace the strong final-depth FT member on 1A2, 2C9 or 3A4.
CheMeleon multiple-depth readouts nearly tie the large 1A2 ensemble but do
not offer a supported gain. The 2C9 and 3A4 selection procedures reject the
new CheMeleon members when rerun. The bootstrap conditions on fitted OOF
predictions and does not account for scenario selection, model refitting or
blind-test shift. No new candidate or submission was prepared.
