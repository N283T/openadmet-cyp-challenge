# D6 gate candidate package: `L2_legacy_x_sqrt_pa_morgan` (alternative)

D6 score: `pct_rank(legacy) * sqrt(gate_x_tdi.p_active)` over the 750 blind compounds; top 49 called.
Legacy = deployed Sep 15 blind vector (`deployed_mean_block`, run `20260924T045344348206Z-7e2db311`); pool and gate from run `20261002T215345688167Z-7920b49a`.
The call vector reproduces run `20261002T221918697661Z-2d8d20f5` `blind_candidates.parquet`.

D6 vs Sep 15 payload: 23 shared, 26 added, 26 removed (boundary gap 0.005919).
D6 vs last sent payload `20261002T092348406577Z-4ea606c1` (175 calls): 28 shared, 21 added, 147 removed.
Mean Morgan gate p_active of called compounds 0.526 (Sep 15: 0.322).
A4: 249 calls, bitwise identical to the last sent payload and to the historical best-scored Sep 15 vector.

## D6 changes versus the Sep 15 payload

| Change | Molecule | Position | Legacy position | Pool rank | p_active |
|---|---|---:|---:|---:|---:|
| added | OCNT-2535689 | 2 | 108 | 0.941 | 0.949 |
| added | OCNT-2535391 | 9 | 91 | 0.869 | 0.641 |
| added | OCNT-2535771 | 10 | 182 | 0.882 | 0.849 |
| added | OCNT-2534867 | 11 | 54 | 0.853 | 0.563 |
| added | OCNT-2535554 | 13 | 51 | 0.963 | 0.529 |
| added | OCNT-2535494 | 19 | 57 | 0.910 | 0.492 |
| added | OCNT-2534798 | 21 | 159 | 0.817 | 0.660 |
| added | OCNT-2535393 | 23 | 138 | 0.942 | 0.592 |
| added | OCNT-2534896 | 24 | 208 | 0.945 | 0.750 |
| added | OCNT-2535292 | 25 | 175 | 0.917 | 0.666 |
| added | OCNT-2535649 | 26 | 120 | 0.722 | 0.554 |
| added | OCNT-2534912 | 28 | 59 | 0.973 | 0.434 |
| added | OCNT-2535489 | 29 | 97 | 0.919 | 0.478 |
| added | OCNT-2534888 | 31 | 131 | 0.975 | 0.507 |
| added | OCNT-2535547 | 34 | 289 | 0.899 | 0.898 |
| added | OCNT-2535624 | 36 | 56 | 0.914 | 0.392 |
| added | OCNT-2535161 | 37 | 77 | 0.943 | 0.413 |
| added | OCNT-2535692 | 40 | 93 | 0.962 | 0.400 |
| added | OCNT-2535422 | 41 | 99 | 0.770 | 0.400 |
| added | OCNT-2535120 | 42 | 178 | 0.951 | 0.501 |
| added | OCNT-2535156 | 43 | 147 | 0.925 | 0.448 |
| added | OCNT-2535160 | 44 | 156 | 0.889 | 0.451 |
| added | OCNT-2535560 | 46 | 114 | 0.757 | 0.388 |
| added | OCNT-2535124 | 47 | 87 | 0.926 | 0.350 |
| added | OCNT-2535757 | 48 | 288 | 0.966 | 0.718 |
| added | OCNT-2534933 | 49 | 62 | 0.967 | 0.322 |
| removed | OCNT-2534824 | 50 | 24 | 0.977 | 0.283 |
| removed | OCNT-2535520 | 51 | 40 | 0.985 | 0.295 |
| removed | OCNT-2535281 | 55 | 43 | 0.978 | 0.267 |
| removed | OCNT-2535143 | 56 | 45 | 0.819 | 0.268 |
| removed | OCNT-2535495 | 64 | 35 | 0.949 | 0.246 |
| removed | OCNT-2535150 | 72 | 33 | 0.946 | 0.233 |
| removed | OCNT-2535151 | 76 | 49 | 0.929 | 0.236 |
| removed | OCNT-2534787 | 88 | 11 | 0.387 | 0.195 |
| removed | OCNT-2535241 | 90 | 42 | 0.865 | 0.205 |
| removed | OCNT-2534781 | 96 | 31 | 0.714 | 0.192 |
| removed | OCNT-2534785 | 100 | 4 | 0.677 | 0.175 |
| removed | OCNT-2535368 | 119 | 47 | 0.719 | 0.169 |
| removed | OCNT-2534866 | 120 | 25 | 0.786 | 0.159 |
| removed | OCNT-2535650 | 126 | 30 | 0.811 | 0.154 |
| removed | OCNT-2534783 | 134 | 13 | 0.354 | 0.140 |
| removed | OCNT-2535010 | 141 | 10 | 0.519 | 0.135 |
| removed | OCNT-2535394 | 148 | 21 | 0.805 | 0.135 |
| removed | OCNT-2534942 | 167 | 48 | 0.990 | 0.132 |
| removed | OCNT-2535064 | 171 | 14 | 0.867 | 0.115 |
| removed | OCNT-2535720 | 183 | 6 | 0.847 | 0.107 |
| removed | OCNT-2535178 | 190 | 27 | 0.582 | 0.110 |
| removed | OCNT-2535126 | 191 | 41 | 0.898 | 0.114 |
| removed | OCNT-2535715 | 252 | 16 | 0.797 | 0.080 |
| removed | OCNT-2535121 | 269 | 34 | 0.498 | 0.080 |
| removed | OCNT-2534962 | 298 | 18 | 0.205 | 0.067 |
| removed | OCNT-2535543 | 367 | 28 | 0.178 | 0.052 |

No model was fitted and nothing was uploaded.
