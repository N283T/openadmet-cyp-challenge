# D6 gate candidate package: `L4_rankmean_legacy_pool` (primary)

D6 score: `0.5 * (pct_rank(legacy) + pct_rank(pool_lgbm.score))` over the 750 blind compounds; top 49 called.
Legacy = deployed Sep 15 blind vector (`deployed_mean_block`, run `20260924T045344348206Z-7e2db311`); pool and gate from run `20261002T215345688167Z-7920b49a`.
The call vector reproduces run `20261002T221918697661Z-2d8d20f5` `blind_candidates.parquet`.

D6 vs Sep 15 payload: 32 shared, 17 added, 17 removed (boundary gap 0.001333).
D6 vs last sent payload `20261002T092348406577Z-4ea606c1` (175 calls): 42 shared, 7 added, 133 removed.
Mean Morgan gate p_active of called compounds 0.368 (Sep 15: 0.322).
A4: 249 calls, bitwise identical to the last sent payload and to the historical best-scored Sep 15 vector.

## D6 changes versus the Sep 15 payload

| Change | Molecule | Position | Legacy position | Pool rank | p_active |
|---|---|---:|---:|---:|---:|
| added | OCNT-2534916 | 18 | 65 | 0.997 | 0.164 |
| added | OCNT-2535554 | 22 | 51 | 0.963 | 0.529 |
| added | OCNT-2534912 | 23 | 59 | 0.973 | 0.434 |
| added | OCNT-2535170 | 25 | 60 | 0.969 | 0.234 |
| added | OCNT-2534933 | 26 | 62 | 0.967 | 0.322 |
| added | OCNT-2535139 | 30 | 64 | 0.961 | 0.275 |
| added | OCNT-2534939 | 31 | 88 | 0.991 | 0.153 |
| added | OCNT-2534870 | 32 | 66 | 0.958 | 0.203 |
| added | OCNT-2535161 | 38 | 77 | 0.943 | 0.413 |
| added | OCNT-2535624 | 40 | 56 | 0.914 | 0.392 |
| added | OCNT-2535692 | 41 | 93 | 0.962 | 0.400 |
| added | OCNT-2535494 | 42 | 57 | 0.910 | 0.492 |
| added | OCNT-2535348 | 44 | 100 | 0.957 | 0.225 |
| added | OCNT-2534938 | 45 | 115 | 0.971 | 0.127 |
| added | OCNT-2535537 | 46 | 53 | 0.883 | 0.261 |
| added | OCNT-2535124 | 47 | 87 | 0.926 | 0.350 |
| added | OCNT-2535509 | 49 | 111 | 0.955 | 0.284 |
| removed | OCNT-2535394 | 57 | 21 | 0.805 | 0.135 |
| removed | OCNT-2535715 | 60 | 16 | 0.797 | 0.080 |
| removed | OCNT-2534924 | 62 | 37 | 0.821 | 0.376 |
| removed | OCNT-2535650 | 63 | 30 | 0.811 | 0.154 |
| removed | OCNT-2535143 | 66 | 45 | 0.819 | 0.268 |
| removed | OCNT-2534866 | 70 | 25 | 0.786 | 0.159 |
| removed | OCNT-2534781 | 95 | 31 | 0.714 | 0.192 |
| removed | OCNT-2534785 | 96 | 4 | 0.677 | 0.175 |
| removed | OCNT-2535368 | 99 | 47 | 0.719 | 0.169 |
| removed | OCNT-2535690 | 108 | 46 | 0.695 | 0.489 |
| removed | OCNT-2535178 | 130 | 27 | 0.582 | 0.110 |
| removed | OCNT-2535010 | 137 | 10 | 0.519 | 0.135 |
| removed | OCNT-2535121 | 159 | 34 | 0.498 | 0.080 |
| removed | OCNT-2534787 | 184 | 11 | 0.387 | 0.195 |
| removed | OCNT-2534783 | 207 | 13 | 0.354 | 0.140 |
| removed | OCNT-2534962 | 280 | 18 | 0.205 | 0.067 |
| removed | OCNT-2535543 | 301 | 28 | 0.178 | 0.052 |

No model was fitted and nothing was uploaded.
