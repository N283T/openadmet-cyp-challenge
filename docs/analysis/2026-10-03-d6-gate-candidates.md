# D6 gate candidates L4 and L2: prepared, not sent

**2026-10-03 (JST).** The reviewer selected two of the blind legacy x gate
vectors inventoried in the [gate concordance analysis](2026-10-03-d6-gate-concordance.md)
(run `20261002T221918697661Z-2d8d20f5`, Part 3) for packaging. Both carry
exactly 49 CYP2D6 positive calls and the exact best-scored 249-call CYP3A4
vector. Neither payload has been uploaded, and no upload is authorized: sending
either prepared file requires explicit user approval for that specific ID.

## Vectors

Ranks are percentile ranks `(rank + 0.5) / 750` from a stable argsort over the
750 blind compounds ordered by `compound_id`; calls are the top 49 under a
stable argsort of the negated score. This reproduces the procedure of run
2d8d20f5, and both packaging runs assert that the resulting call vector equals
that run's `blind_candidates.parquet` vector for the variant.

- Legacy score: the deployed Sep 15 blind vector (`deployed_mean_block` of the
  `blind-scores.parquet` input of run `20260924T045344348206Z-7e2db311`). The
  packaging run asserts that its top 49 reproduce the 49 submitted Sep 15 calls.
- Pool and gate scores: the frozen blind scores of run
  `20261002T215345688167Z-7920b49a` (`pool_lgbm.score`, `gate_x_tdi.p_active`).

| Role | Variant | D6 score | Packaging run | Candidate | Prepared submission |
|---|---|---|---|---|---|
| Primary | `L4_rankmean_legacy_pool` | `0.5 * (rank(legacy) + rank(pool_lgbm))` | `20261002T223203072149Z-662a599c` | `20261002T223224277051Z-1168833a` | `20261002T223229553466Z-a08301d2` |
| Alternative | `L2_legacy_x_sqrt_pa_morgan` | `rank(legacy) * sqrt(p_active)` | `20261002T223215615405Z-0def2858` | `20261002T223224511229Z-f62f5093` | `20261002T223229749698Z-16026b6c` |

Call-boundary gaps are .00133 (L4) and .00592 (L2); neither selection is tied.

CYP3A4 is copied unchanged from the last sent payload
`20261002T092348406577Z-4ea606c1`. Its 249 calls are bitwise identical to the
CYP3A4 column of the Sep 15 payload as re-prepared in
`20260924T050705396203Z-10e1d3ad`, which shares the SHA-256 of the originally
sent `20260915T130004513006Z-75b25070` file and is the historical best-scored A4
vector. Both packaging runs and an independent read of the prepared files
confirm this identity.

## Counts

| Prepared ID | D6 calls | vs Sep 15 (49): shared / added / removed | vs last sent 4ea606c1 (175): shared / added / removed | A4 calls | A4 changed vs last sent | A4 changed vs Sep 15 |
|---|---:|---|---|---:|---:|---:|
| `20261002T223229553466Z-a08301d2` (L4) | 49 | 32 / 17 / 17 | 42 / 7 / 133 | 249 | 0 | 0 |
| `20261002T223229749698Z-16026b6c` (L2) | 49 | 23 / 26 / 26 | 28 / 21 / 147 | 249 | 0 | 0 |

## Dry-run validation

`pixi run cyp submit ID --dry-run` succeeded for both prepared files: 750 rows,
the four expected columns, valid identifiers and a complete private profile.
No live submit command was issued.

| Prepared ID | SHA-256 | Target changes vs last sent | Next eligible |
|---|---|---|---|
| `20261002T223229553466Z-a08301d2` | `f6e2ddb44f86ee11aea1877ff07dc79a737f54cb7f9793e63b924af7893af393` | D6 140, A4 0 | 2026-10-03T09:08:22+09:00 |
| `20261002T223229749698Z-16026b6c` | `48c57605cca07396e368abacebc66b38a205cdec7a5158d0f7b0ec10879c3c5d` | D6 168, A4 0 | 2026-10-03T09:08:22+09:00 |

The classification track was last sent on 2026-10-02 at 21:08:22 JST, so the
12-hour cooldown ends at the time above. One upload of either payload consumes
that opportunity for the track.

## Changed-call diagnostic

Run `20261002T223515148996Z-ccb05cf3` is read-only over the two packaging-run
score tables and the two prepared files; it also asserts that the recorded
p_active and pool scores equal the `gate_x_tdi` and `pool_lgbm` rows of
`test_scores.parquet` in run 7920b49a. Calls are grouped against the Sep 15
payload as retained (called by both), removed (Sep 15 only) and added (candidate
only).

The neighbour rate is the similarity-weighted mean label of the 10 most similar
labelled D6 training molecules, where labelled means carrying a `CYP2D6_is_TDI`
measurement (1,497 molecules, 324 positive). Similarity is Tanimoto on Morgan
radius-2 2,048-bit fingerprints of the standardized SMILES of the canonical
`compounds.parquet`. Positions are 1-based with 1 the highest score; combined
refers to the candidate's own score.

### L4_rankmean_legacy_pool group means

| Group | n | Legacy rank | Pool rank | Combined rank | Legacy pos. | Pool pos. | Combined pos. | p_active | Top-10 weighted pos. rate | Nearest pos. sim. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| retained | 32 | 0.968 | 0.948 | 0.975 | 24.4 | 39.7 | 19.6 | 0.399 | 0.466 | 0.460 |
| removed | 17 | 0.966 | 0.610 | 0.830 | 26.1 | 293.1 | 127.9 | 0.175 | 0.348 | 0.392 |
| added | 17 | 0.902 | 0.953 | 0.954 | 74.4 | 35.8 | 35.2 | 0.309 | 0.393 | 0.463 |

### L4_rankmean_legacy_pool changed calls

| Change | Molecule | Legacy rank (pos.) | Pool rank (pos.) | Combined rank (pos.) | p_active | Top-10 weighted pos. rate |
|---|---|---:|---:|---:|---:|---:|
| removed | OCNT-2535394 | 0.973 (21) | 0.805 (147) | 0.923 (57) | 0.135 | 0.441 |
| removed | OCNT-2535715 | 0.979 (16) | 0.797 (153) | 0.919 (60) | 0.080 | 0.426 |
| removed | OCNT-2534924 | 0.951 (37) | 0.821 (135) | 0.917 (62) | 0.376 | 0.079 |
| removed | OCNT-2535650 | 0.961 (30) | 0.811 (142) | 0.918 (63) | 0.154 | 0.388 |
| removed | OCNT-2535143 | 0.941 (45) | 0.819 (136) | 0.913 (66) | 0.268 | 0.525 |
| removed | OCNT-2534866 | 0.967 (25) | 0.786 (161) | 0.907 (70) | 0.159 | 0.335 |
| removed | OCNT-2534781 | 0.959 (31) | 0.714 (215) | 0.875 (95) | 0.192 | 0.342 |
| removed | OCNT-2534785 | 0.995 (4) | 0.677 (243) | 0.873 (96) | 0.175 | 0.492 |
| removed | OCNT-2535368 | 0.938 (47) | 0.719 (211) | 0.869 (99) | 0.169 | 0.326 |
| removed | OCNT-2535690 | 0.939 (46) | 0.695 (229) | 0.857 (108) | 0.489 | 0.128 |
| removed | OCNT-2535178 | 0.965 (27) | 0.582 (314) | 0.827 (130) | 0.110 | 0.428 |
| removed | OCNT-2535010 | 0.987 (10) | 0.519 (361) | 0.818 (137) | 0.135 | 0.353 |
| removed | OCNT-2535121 | 0.955 (34) | 0.498 (377) | 0.789 (159) | 0.080 | 0.588 |
| removed | OCNT-2534787 | 0.986 (11) | 0.387 (460) | 0.755 (184) | 0.195 | 0.356 |
| removed | OCNT-2534783 | 0.983 (13) | 0.354 (485) | 0.725 (207) | 0.140 | 0.349 |
| removed | OCNT-2534962 | 0.977 (18) | 0.205 (597) | 0.627 (280) | 0.067 | 0.274 |
| removed | OCNT-2535543 | 0.963 (28) | 0.178 (617) | 0.599 (301) | 0.052 | 0.093 |
| added | OCNT-2534916 | 0.915 (65) | 0.997 (3) | 0.977 (18) | 0.164 | 0.415 |
| added | OCNT-2535554 | 0.933 (51) | 0.963 (28) | 0.971 (22) | 0.529 | 0.651 |
| added | OCNT-2534912 | 0.922 (59) | 0.973 (21) | 0.969 (23) | 0.434 | 0.533 |
| added | OCNT-2535170 | 0.921 (60) | 0.969 (24) | 0.967 (25) | 0.234 | 0.268 |
| added | OCNT-2534933 | 0.918 (62) | 0.967 (25) | 0.966 (26) | 0.322 | 0.302 |
| added | OCNT-2535139 | 0.914 (64) | 0.961 (30) | 0.959 (30) | 0.275 | 0.504 |
| added | OCNT-2534939 | 0.883 (88) | 0.991 (7) | 0.961 (31) | 0.153 | 0.461 |
| added | OCNT-2534870 | 0.913 (66) | 0.958 (32) | 0.958 (32) | 0.203 | 0.468 |
| added | OCNT-2535161 | 0.898 (77) | 0.943 (43) | 0.950 (38) | 0.413 | 0.506 |
| added | OCNT-2535624 | 0.926 (56) | 0.914 (65) | 0.949 (40) | 0.392 | 0.294 |
| added | OCNT-2535692 | 0.877 (93) | 0.962 (29) | 0.946 (41) | 0.400 | 0.256 |
| added | OCNT-2535494 | 0.925 (57) | 0.910 (68) | 0.945 (42) | 0.492 | 0.549 |
| added | OCNT-2535348 | 0.867 (100) | 0.957 (33) | 0.942 (44) | 0.225 | 0.487 |
| added | OCNT-2534938 | 0.847 (115) | 0.971 (22) | 0.941 (45) | 0.127 | 0.416 |
| added | OCNT-2535537 | 0.930 (53) | 0.883 (88) | 0.939 (46) | 0.261 | 0.238 |
| added | OCNT-2535124 | 0.885 (87) | 0.926 (56) | 0.938 (47) | 0.350 | 0.240 |
| added | OCNT-2535509 | 0.853 (111) | 0.955 (34) | 0.935 (49) | 0.284 | 0.096 |

### L2_legacy_x_sqrt_pa_morgan group means

| Group | n | Legacy rank | Pool rank | Combined rank | Legacy pos. | Pool pos. | Combined pos. | p_active | Top-10 weighted pos. rate | Nearest pos. sim. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| retained | 23 | 0.972 | 0.940 | 0.975 | 21.3 | 45.7 | 19.4 | 0.496 | 0.430 | 0.469 |
| removed | 26 | 0.963 | 0.734 | 0.812 | 28.3 | 200.1 | 141.3 | 0.167 | 0.421 | 0.407 |
| added | 26 | 0.833 | 0.904 | 0.961 | 126.0 | 72.7 | 29.9 | 0.552 | 0.317 | 0.387 |

### L2_legacy_x_sqrt_pa_morgan changed calls

| Change | Molecule | Legacy rank (pos.) | Pool rank (pos.) | Combined rank (pos.) | p_active | Top-10 weighted pos. rate |
|---|---|---:|---:|---:|---:|---:|
| removed | OCNT-2534824 | 0.969 (24) | 0.977 (18) | 0.934 (50) | 0.283 | 0.269 |
| removed | OCNT-2535520 | 0.947 (40) | 0.985 (12) | 0.933 (51) | 0.295 | 0.605 |
| removed | OCNT-2535281 | 0.943 (43) | 0.978 (17) | 0.927 (55) | 0.267 | 0.552 |
| removed | OCNT-2535143 | 0.941 (45) | 0.819 (136) | 0.926 (56) | 0.268 | 0.525 |
| removed | OCNT-2535495 | 0.954 (35) | 0.949 (39) | 0.915 (64) | 0.246 | 0.445 |
| removed | OCNT-2535150 | 0.957 (33) | 0.946 (41) | 0.905 (72) | 0.233 | 0.632 |
| removed | OCNT-2535151 | 0.935 (49) | 0.929 (54) | 0.899 (76) | 0.236 | 0.453 |
| removed | OCNT-2534787 | 0.986 (11) | 0.387 (460) | 0.883 (88) | 0.195 | 0.356 |
| removed | OCNT-2535241 | 0.945 (42) | 0.865 (102) | 0.881 (90) | 0.205 | 0.272 |
| removed | OCNT-2534781 | 0.959 (31) | 0.714 (215) | 0.873 (96) | 0.192 | 0.342 |
| removed | OCNT-2534785 | 0.995 (4) | 0.677 (243) | 0.867 (100) | 0.175 | 0.492 |
| removed | OCNT-2535368 | 0.938 (47) | 0.719 (211) | 0.842 (119) | 0.169 | 0.326 |
| removed | OCNT-2534866 | 0.967 (25) | 0.786 (161) | 0.841 (120) | 0.159 | 0.335 |
| removed | OCNT-2535650 | 0.961 (30) | 0.811 (142) | 0.833 (126) | 0.154 | 0.388 |
| removed | OCNT-2534783 | 0.983 (13) | 0.354 (485) | 0.822 (134) | 0.140 | 0.349 |
| removed | OCNT-2535010 | 0.987 (10) | 0.519 (361) | 0.813 (141) | 0.135 | 0.353 |
| removed | OCNT-2535394 | 0.973 (21) | 0.805 (147) | 0.803 (148) | 0.135 | 0.441 |
| removed | OCNT-2534942 | 0.937 (48) | 0.990 (8) | 0.778 (167) | 0.132 | 0.472 |
| removed | OCNT-2535064 | 0.982 (14) | 0.867 (100) | 0.773 (171) | 0.115 | 0.391 |
| removed | OCNT-2535720 | 0.993 (6) | 0.847 (115) | 0.757 (183) | 0.107 | 0.622 |
| removed | OCNT-2535178 | 0.965 (27) | 0.582 (314) | 0.747 (190) | 0.110 | 0.428 |
| removed | OCNT-2535126 | 0.946 (41) | 0.898 (77) | 0.746 (191) | 0.114 | 0.510 |
| removed | OCNT-2535715 | 0.979 (16) | 0.797 (153) | 0.665 (252) | 0.080 | 0.426 |
| removed | OCNT-2535121 | 0.955 (34) | 0.498 (377) | 0.642 (269) | 0.080 | 0.588 |
| removed | OCNT-2534962 | 0.977 (18) | 0.205 (597) | 0.603 (298) | 0.067 | 0.274 |
| removed | OCNT-2535543 | 0.963 (28) | 0.178 (617) | 0.511 (367) | 0.052 | 0.093 |
| added | OCNT-2535689 | 0.857 (108) | 0.941 (45) | 0.998 (2) | 0.949 | 0.451 |
| added | OCNT-2535391 | 0.879 (91) | 0.869 (99) | 0.989 (9) | 0.641 | 0.218 |
| added | OCNT-2535771 | 0.758 (182) | 0.882 (89) | 0.987 (10) | 0.849 | 0.290 |
| added | OCNT-2534867 | 0.929 (54) | 0.853 (111) | 0.986 (11) | 0.563 | 0.265 |
| added | OCNT-2535554 | 0.933 (51) | 0.963 (28) | 0.983 (13) | 0.529 | 0.651 |
| added | OCNT-2535494 | 0.925 (57) | 0.910 (68) | 0.975 (19) | 0.492 | 0.549 |
| added | OCNT-2534798 | 0.789 (159) | 0.817 (138) | 0.973 (21) | 0.660 | 0.437 |
| added | OCNT-2535393 | 0.817 (138) | 0.942 (44) | 0.970 (23) | 0.592 | 0.324 |
| added | OCNT-2534896 | 0.723 (208) | 0.945 (42) | 0.969 (24) | 0.750 | 0.295 |
| added | OCNT-2535292 | 0.767 (175) | 0.917 (63) | 0.967 (25) | 0.666 | 0.294 |
| added | OCNT-2535649 | 0.841 (120) | 0.722 (209) | 0.966 (26) | 0.554 | 0.290 |
| added | OCNT-2534912 | 0.922 (59) | 0.973 (21) | 0.963 (28) | 0.434 | 0.533 |
| added | OCNT-2535489 | 0.871 (97) | 0.919 (61) | 0.962 (29) | 0.478 | 0.463 |
| added | OCNT-2534888 | 0.826 (131) | 0.975 (19) | 0.959 (31) | 0.507 | 0.426 |
| added | OCNT-2535547 | 0.615 (289) | 0.899 (76) | 0.955 (34) | 0.898 | 0.132 |
| added | OCNT-2535624 | 0.926 (56) | 0.914 (65) | 0.953 (36) | 0.392 | 0.294 |
| added | OCNT-2535161 | 0.898 (77) | 0.943 (43) | 0.951 (37) | 0.413 | 0.506 |
| added | OCNT-2535692 | 0.877 (93) | 0.962 (29) | 0.947 (40) | 0.400 | 0.256 |
| added | OCNT-2535422 | 0.869 (99) | 0.770 (173) | 0.946 (41) | 0.400 | 0.096 |
| added | OCNT-2535120 | 0.763 (178) | 0.951 (37) | 0.945 (42) | 0.501 | 0.314 |
| added | OCNT-2535156 | 0.805 (147) | 0.925 (57) | 0.943 (43) | 0.448 | 0.187 |
| added | OCNT-2535160 | 0.793 (156) | 0.889 (84) | 0.942 (44) | 0.451 | 0.092 |
| added | OCNT-2535560 | 0.849 (114) | 0.757 (183) | 0.939 (46) | 0.388 | 0.091 |
| added | OCNT-2535124 | 0.885 (87) | 0.926 (56) | 0.938 (47) | 0.350 | 0.240 |
| added | OCNT-2535757 | 0.617 (288) | 0.966 (26) | 0.937 (48) | 0.718 | 0.233 |
| added | OCNT-2534933 | 0.918 (62) | 0.967 (25) | 0.935 (49) | 0.322 | 0.302 |

### Interpretation

These group means describe the composition of the swaps; they are not estimates
of the hidden labels and carry no validated performance claim.

L4 drops Sep 15 calls that `pool_lgbm` ranks far down (mean pool position 293,
mean p_active .175) and promotes compounds both scores rank highly (mean legacy
position 74, mean pool position 36). Its added set has somewhat more
labelled-neighbour support than its removed set (top-10 weighted positive rate
.393 against .348; mean nearest-positive similarity .463 against .392), and the
32 retained calls remain the best-supported group (.466).

L2 pushes the call set toward high gate probability (added mean p_active .552
against removed .167) at the cost of legacy ranking (added mean legacy position
126) and of labelled-neighbour support: its added set has the lowest neighbour
positive rate of any group here (.317 against .421 removed and .430 retained),
and several additions sit near .09-.13. Four additions come from legacy
positions beyond 180. This is the same label-set tension recorded in RL0213 for
the product forms, now visible in the blind composition.

Neither diagnostic establishes either vector as an improvement. The gate's
directional agreement with historical Live movements (RL0213) is the only
out-of-sample evidence behind these swaps, it rests on seven overlapping
transitions against one hidden set, and no variant has a validated
non-conflicting OOF gain.

## Status

No upload has occurred through `cyp submit`, the Space UI or API, or any other
route; both prepared files remain local with `uploaded: false` recorded in their
packaging summaries. Preparation and dry-running do not authorize an upload.
Sending either payload requires explicit user approval for that specific
prepared ID, immediately beforehand, and would consume the 12-hour
classification opportunity. User approval is pending, and the two payloads are
mutually exclusive for the next opportunity.

Compact outputs, provenance and the dry-run records are promoted to
`artifacts/tracked/cyp-d6-gate-candidate-20261002T223203072149Z-662a599c/`,
`artifacts/tracked/cyp-d6-gate-candidate-20261002T223215615405Z-0def2858/` and
`artifacts/tracked/cyp-d6-gate-candidate-diagnostic-20261002T223515148996Z-ccb05cf3/`.
