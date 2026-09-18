# Regression foundation-readout v8 candidate

## Screen design

Run `3101f488` evaluates TabPFN 3.5-fast and TabICL readouts over twelve
target-specific feature combinations selected from the preceding LightGBM and
dual-split screens. High-dimensional pool selection is repeated inside every
outer training fold with the same fixed LightGBM selector; held-out rows do not
select features. TabICL receives outer-training median imputation and scaling,
while TabPFN receives outer-training median imputation and its internal
preprocessing.

The useful primary-split results are:

| Target / features | LGBM | TabPFN fast | TabICL |
|---|---:|---:|---:|
| 1A2 all pool top512 | .693368 | **.690505** | .703441 |
| 2C9 dense pool top512 | .496386 | **.488792** | .489547 |
| 3A4 dense pool top512 | .423791 | **.423025** | .423723 |

CYP2D6 foundation readouts range from `.843974` to `.874744` and are rejected.
Full TabPFN 3.5 run `1ff1c74c` is worse than fast on all four promoted
target/view pairs, so fast/full duplicates are not carried together.

Scaffold confirmation `79d216a4` gives `.702839` for 1A2 all-pool TabPFN,
`.495002/.499007` for 2C9 dense-pool TabPFN/TabICL, and
`.407138/.406291` for 3A4. The 3A4 foundation readouts improve its scaffold
LightGBM `.413217`; 2C9 remains in the leading group but below scaffold LightGBM
`.491042`; the 1A2 single-model improvement does not transfer to scaffold.

Same-feature LGBM versus foundation Spearman correlations are `.9484` for 1A2,
`.9707/.9718` for 2C9 and `.9793/.9788` for 3A4. Fast versus TabICL is much
more redundant at `.9797/.9916/.9951`. Thus no target retains both new
foundation backends merely as duplicates.

## Bank review

Run `46314b50` compares v7, LGBM swaps, individual foundation additions and
joint LGBM/foundation additions using the unchanged flat-bank methods.

| Target | v7 | Selected v8 | Method | Change |
|---|---:|---:|---|---:|
| 1A2 | .674719 | **.670833** | ten-member equal mean | -.003886 |
| 2C9 | .445372 | **.444776** | seven-member interval-simplex | -.000596 |
| 2D6 | .814706 | **.814706** | thirteen-member Caruana | unchanged |
| 3A4 | .387486 | **.386084** | eleven-member equal mean | -.001402 |

The 1A2 bank keeps both prior LGBM pools and adds all-pool TabPFN because the
three predictions are not near duplicates. The 2C9 bank adds dense-pool TabICL;
adding both dense LGBM and a foundation readout does not improve it. The 3A4
bank benefits from retaining dense top256, adding dense top512 LightGBM, and
adding the fast foundation readout. Equal means win cleanly for 1A2 and 3A4,
avoiding post-hoc constant weights.

## Full fit and candidate

Run `a105b341` full-fits eleven newly required atomic members using all available
endpoint labels and predicts the 750 blinded compounds. Run `21ad1846` supplies
the missing full CheMeleon TabPFN component for 2C9. Existing atomic test
predictions are reused from the last submitted candidate only when the atomic
recipe is unchanged.

Packaging run `b97ec048` creates both raw and placement-aligned candidates. The
registered candidate is `1f52e80c`, using the previous submitted target mean and
standard deviation without changing ranks. All four targets change. Relative to
the previous submission, placed mean absolute changes are `.09957/.06687/.04704/
.07270`, with rank correlations `.99352/.99547/.99193/.99562` for
1A2/2C9/D6/3A4. No test label is used. CYP2D6 remains the risky component because
its previous OOF gain failed to transfer; this candidate is prepared for review,
not uploaded.

