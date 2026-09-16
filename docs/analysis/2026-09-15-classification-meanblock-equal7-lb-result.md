# D6 mean-block plus A4 equal-seven leaderboard result

## Result

The September 15 13:00 UTC classification submission is reflected on all three
leaderboards. It raises macro MCC from 0.36409 (rank 13) to 0.41927 (rank 3).

| Endpoint | Previous MCC | New MCC | Change | Previous rank | New rank |
|---|---:|---:|---:|---:|---:|
| CYP2D6 TDI | 0.27857 | 0.30864 | +0.03006 | 20 | 15 |
| CYP3A4 TDI | 0.44962 | 0.52991 | +0.08029 | 16 | 1 |
| Macro mean | 0.36409 | 0.41927 | +0.05518 | 13 | 3 |

Both endpoint call rates were held fixed. D6 retained 40 of 49 previous positive
calls and swapped 9 each way. A4 retained 231 of 249 and swapped 18 each way.
At those unchanged rates, precision and recall both increased: D6 precision/recall
rose 0.31967/0.34685 to 0.34121/0.38349; A4 rose 0.57560/0.67354 to
0.62845/0.72921. The gain therefore reflects better selected calls, not a
prevalence or threshold-rate change.

## Attribution and uncertainty

A4 contributes 72.8% of the macro gain and D6 27.2%. This supports both endpoint
changes as a package, but it does not identify individual model contributions.
The A4 submission simultaneously changes the old 80/20 weighting and adds the raw
Jeremy views to a seven-member equal mean. One blind result cannot attribute the
gain specifically to raw embeddings.

The A4 grouped-CV target-MCC improvement was only +0.01074 with a 95% interval of
[-0.00787, 0.02930], whereas its leaderboard gain is +0.08029. The D6 PCA-to-mean
block comparison had identical target-MCC calls and only +0.00018 AUROC in CV.
The larger leaderboard movements may reflect favorable evaluation-set alignment,
ensemble changes outside the isolated raw/PCA contrast, or selection luck. Do not
treat their ratio as a transferable CV-to-leaderboard scaling rule.

The A4 lead over rank 2 is 0.01446, smaller than the reported marginal bootstrap
standard deviations. Those deviations are not paired uncertainty for rank
differences, so they neither prove nor disprove a stable lead. Macro rank 3 is
only 0.01014 above rank 4 and 0.02148 below rank 2.

## Decision

Freeze this exact A4 vector as the live reference. Do not immediately perturb it
or claim that raw embeddings alone caused the gain. With A4 held at 0.52991, D6
would need MCC 0.35200 to match the current rank-1 macro score, an improvement of
0.04337 from 0.30864. D6 is therefore the highest-value next endpoint, while A4
ablation remains a scientific analysis rather than the next live submission.

Recorded analysis run: `20260915T131233106223Z-d6e20fb5`.
