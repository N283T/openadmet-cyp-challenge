# Initial descriptor correlation findings — 2026-09-09

Analysis run: `20260909T133949502832Z-d7b6d4e2`.

Regression shows descriptive signal: MolLogP has Spearman rho 0.600 for CYP3A4
(n=2,335), and 0.485 for CYP2C9 (n=1,285). CYP2D6 instead has nBase as the top
single descriptor (rho 0.305); its MolLogP correlation is only 0.034.
For CYP1A2, ATS3are ranks first (rho -0.367); the more interpretable FractionCSP3
also correlates negatively (rho -0.309).

Binary TDI correlations are much weaker: the maximum eligible absolute rho is
0.110 for CYP2D6 (nS; n=1,497), and 0.179 for CYP3A4 (QED; n=3,584).
These binary-target correlations are not MCC scores or measures of test accuracy.
Positive label counts are 324/1,497 and 764/3,584, respectively.

Restricting measured pIC50 to >= 4 reduces MolLogP rho from 0.600 to 0.374 for
CYP3A4 (n=1,391), and 0.485 to 0.356 for CYP2C9 (n=1,026). This sensitivity
analysis conditions on the outcome and changes its range; it is not a proposed
training filter. Raw full-range correlations partly reflect separation of weak
inhibitors from more active compounds. The point-estimate correlations do not
reproduce the competition's interval-aware ST-RAE metric.

Interpretation: useful univariate associations exist in regression, while these
simple descriptors do not strongly separate TDI classes individually. Many top
features describe similar properties, so they should not be counted as independent
signals. No evidence here establishes causality, generalization to the analog test
set, or that nonlinear/multivariate TDI models will fail.

Validation: standardized molecule/endpoint uniqueness and join cardinality checked;
no test labels used; pairwise missingness recorded; top Pearson coefficients
independently reproduced with NumPy. Bootstrap intervals and BH-adjusted p-values
are exploratory because descriptors were selected on the same data and related
molecules are not independent. No split, model, imputation, or feature selection.

Full correlations, common descriptors, sample counts and conditional bootstrap
intervals are under the analysis run's `outputs/` directory. Source code:
`scripts/descriptor_correlations.py`; config: `configs/analysis/correlations.json`.

---

# Descriptor–target correlation screening

Source dataset: `20260909T132734119233Z-a7063ccd`; feature set: `20260909T132918829974Z-bc81af58`.
Measured training labels only, one row per standardized compound and endpoint.
Pairwise complete cases, no imputation, models, splits, or feature removal.
Rankings require n >= 100 and >= 80% descriptor coverage. All tested pairs are in correlations.csv.
BH q-values span all valid endpoint/descriptor pairs; molecular dependence limits inference.
This is exploratory selection on the same data, not out-of-sample performance.
Top-feature bootstrap intervals are conditional on selection and assume independent molecules.

## CYP1A2_pIC50_direct_inhibition

| Descriptor | n | Coverage | Pearson | Spearman |
|---|---:|---:|---:|---:|
| mordred:ATS3are | 1412 | 100.0% | -0.341 | -0.367 |
| mordred:ATS3pe | 1412 | 100.0% | -0.338 | -0.365 |
| mordred:ATS3se | 1412 | 100.0% | -0.336 | -0.364 |
| mordred:ATS3i | 1412 | 100.0% | -0.333 | -0.364 |
| mordred:nBondsS | 1412 | 100.0% | -0.329 | -0.363 |
| mordred:ATS2i | 1412 | 100.0% | -0.324 | -0.363 |
| mordred:ATS4are | 1412 | 100.0% | -0.335 | -0.362 |
| mordred:ATS4pe | 1412 | 100.0% | -0.333 | -0.360 |
| mordred:ATS4se | 1412 | 100.0% | -0.332 | -0.359 |
| mordred:ATS2se | 1412 | 100.0% | -0.319 | -0.355 |

## CYP2C9_pIC50_direct_inhibition

| Descriptor | n | Coverage | Pearson | Spearman |
|---|---:|---:|---:|---:|
| rdkit:MolLogP | 1285 | 100.0% | 0.469 | 0.485 |
| mordred:SLogP | 1285 | 100.0% | 0.469 | 0.485 |
| mordred:FilterItLogS | 1285 | 100.0% | -0.463 | -0.442 |
| mordred:nAromAtom | 1285 | 100.0% | 0.390 | 0.393 |
| mordred:nAromBond | 1285 | 100.0% | 0.383 | 0.389 |
| mordred:nBondsA | 1285 | 100.0% | 0.383 | 0.389 |
| mordred:ETA_beta_ns | 1285 | 100.0% | 0.384 | 0.382 |
| rdkit:SMR_VSA7 | 1285 | 100.0% | 0.380 | 0.374 |
| mordred:SMR_VSA7 | 1285 | 100.0% | 0.380 | 0.374 |
| mordred:C2SP2 | 1285 | 100.0% | 0.362 | 0.369 |

## CYP2D6_pIC50_direct_inhibition

| Descriptor | n | Coverage | Pearson | Spearman |
|---|---:|---:|---:|---:|
| mordred:nBase | 1493 | 100.0% | 0.265 | 0.305 |
| mordred:AATS4s | 1493 | 100.0% | -0.232 | -0.236 |
| mordred:AATS4se | 1493 | 100.0% | -0.231 | -0.235 |
| mordred:GATS2pe | 1493 | 100.0% | 0.242 | 0.231 |
| mordred:AATS4are | 1493 | 100.0% | -0.228 | -0.229 |
| mordred:AATS5pe | 1493 | 100.0% | -0.220 | -0.229 |
| mordred:AATS5se | 1493 | 100.0% | -0.216 | -0.229 |
| mordred:AATS5are | 1493 | 100.0% | -0.216 | -0.228 |
| mordred:GATS2are | 1493 | 100.0% | 0.237 | 0.226 |
| mordred:AATS4pe | 1493 | 100.0% | -0.229 | -0.226 |

## CYP3A4_pIC50_direct_inhibition

| Descriptor | n | Coverage | Pearson | Spearman |
|---|---:|---:|---:|---:|
| rdkit:MolLogP | 2335 | 100.0% | 0.571 | 0.600 |
| mordred:SLogP | 2335 | 100.0% | 0.571 | 0.600 |
| mordred:FilterItLogS | 2335 | 100.0% | -0.526 | -0.550 |
| mordred:MINdO | 2135 | 91.4% | 0.425 | 0.394 |
| mordred:nAromAtom | 2335 | 100.0% | 0.362 | 0.393 |
| mordred:ATS0p | 2335 | 100.0% | 0.360 | 0.388 |
| mordred:nAromBond | 2335 | 100.0% | 0.355 | 0.388 |
| mordred:nBondsA | 2335 | 100.0% | 0.355 | 0.388 |
| mordred:C2SP2 | 2335 | 100.0% | 0.375 | 0.384 |
| mordred:ETA_beta_ns | 2335 | 100.0% | 0.357 | 0.380 |

## CYP2D6_is_TDI

| Descriptor | n | Coverage | Pearson | Spearman |
|---|---:|---:|---:|---:|
| mordred:nS | 1497 | 100.0% | 0.107 | 0.110 |
| mordred:BCUTi-1l | 1497 | 100.0% | -0.109 | -0.096 |
| mordred:n9FHRing | 1497 | 100.0% | 0.094 | 0.091 |
| mordred:SaaS | 1497 | 100.0% | 0.090 | 0.090 |
| mordred:NaaS | 1497 | 100.0% | 0.088 | 0.087 |
| rdkit:PEOE_VSA12 | 1497 | 100.0% | 0.070 | 0.087 |
| mordred:PEOE_VSA12 | 1497 | 100.0% | 0.070 | 0.087 |
| mordred:GATS7c | 1497 | 100.0% | -0.084 | -0.084 |
| mordred:AATSC0p | 1497 | 100.0% | 0.054 | 0.082 |
| rdkit:fr_furan | 1497 | 100.0% | 0.082 | 0.082 |

## CYP3A4_is_TDI

| Descriptor | n | Coverage | Pearson | Spearman |
|---|---:|---:|---:|---:|
| rdkit:qed | 3584 | 100.0% | 0.185 | 0.179 |
| mordred:GhoseFilter | 3584 | 100.0% | 0.160 | 0.160 |
| mordred:Xch-5d | 3584 | 100.0% | 0.151 | 0.136 |
| mordred:Kier2 | 3584 | 100.0% | -0.142 | -0.128 |
| mordred:Xch-4dv | 3584 | 100.0% | 0.137 | 0.128 |
| mordred:Xch-4d | 3584 | 100.0% | 0.137 | 0.128 |
| rdkit:Phi | 3584 | 100.0% | -0.134 | -0.126 |
| mordred:VR2_A | 3584 | 100.0% | -0.070 | -0.123 |
| mordred:Kier3 | 3584 | 100.0% | -0.119 | -0.122 |
| mordred:AXp-1d | 3584 | 100.0% | -0.123 | -0.121 |
