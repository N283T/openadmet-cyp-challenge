# Frozen CYP/HEM trunk: labelled-only TDI integration diagnostics

Diagnostic source run: `20261001T133225940252Z-502d403a`. No new model fits or trunk extractions.

The user explicitly deferred full-pool extraction because of cost, effort and storage. This analysis looks for signals worth completing the operating pool later. It cannot set a deployment rate or justify a submitted call vector. Regression is unchanged.

## Cohort and controls

| target | labelled_rows | total_pool_rows | missing_rows | reconstruction_max_diff |
| --- | --- | --- | --- | --- |
| d6 | 1497 | 6146 | 4649 | 0.000000 |
| a4 | 3584 | 6146 | 2562 | 0.000000 |

Each trunk recipe covers exactly the incumbent labelled molecule IDs and UMAP folds. Measurement labels and saved fold assignments are checked explicitly. The exact saved incumbent composition is reproduced. Missing classification labels are never negatives.

D6 preserves its saved core-member choices, PCA/raw TabPFN block, rank transform and Morgan logistic/SVM/complement-NB tail stages. A4 preserves its seven rank members: historical stacked TabPFN with Mordred, pretrained molecular 600-dimensional embeddings with PCA256/raw600 and TabICL/TabPFN, external-CYP classification embeddings plus Mordred with LightGBM, and Morgan plus PXR/ADMET auxiliary predictions with LightGBM.

The no-trunk control uses the same saved original scores and tests adding an existing original atom as a core/equal-mean member; it diagnoses weight/composition changes. It is not an unrestricted reselection of every historical model bank. Molecular-only CheMeleon2048 + Mordred1613 readouts are separately included as ancillary placements.

## Conditional selection procedures

Single addition/replacement is chosen on the other four saved OOF folds, separately by AP and target-prior MCC at the historical diagnostic rate. The joint policy selects both single placement and rate on those four folds. All include keeping the incumbent.

Multi-block selection uses inner held folds among those four to select block sizes 1/2/3/5 and either one averaged member, multiple individual additions, or 10/25/50% mixing of incumbent ranks with trunk-block ranks. Distinct readouts are greedily selected on inner training folds, then reselected on all four folds for the outer score. Greedy replacement preserves existing member positions and stops when training-fold quality no longer improves: D6 replaces only its PCA/raw block; A4 can replace up to three distinct members with distinct readouts.

All methods are conditional saved-OOF analyses, not fully nested base-model refits. Other-fold OOF base models may have used evaluated-fold labels in their original training. Multiple trunks here are readouts/views of the same frozen extraction, not independently computed Boltz trunks or receptor conformations.

## Selected policies

| target | scenario | AP | AUROC | MCC_target_prior | added | removed | positive_added | positive_removed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d6 | current | 0.296860 | 0.580764 | 0.110109 | 0 | 0 | 0 | 0 |
| d6 | single_selected_AP | 0.306124 | 0.591360 | 0.099878 | 9 | 9 | 1 | 3 |
| d6 | no_trunk_selected_AP | 0.293968 | 0.581306 | 0.104965 | 12 | 12 | 3 | 4 |
| d6 | multi_block_selected_AP | 0.318808 | 0.611664 | 0.104965 | 45 | 45 | 17 | 18 |
| d6 | multi_replace_selected_AP | 0.299876 | 0.585624 | 0.099878 | 9 | 9 | 1 | 3 |
| d6 | single_selected_MCC | 0.300191 | 0.590498 | 0.099878 | 12 | 12 | 3 | 5 |
| d6 | no_trunk_selected_MCC | 0.295221 | 0.577591 | 0.110109 | 7 | 7 | 1 | 1 |
| d6 | multi_block_selected_MCC | 0.302266 | 0.585844 | 0.104965 | 19 | 19 | 6 | 7 |
| d6 | multi_replace_selected_MCC | 0.298569 | 0.579246 | 0.099878 | 2 | 2 | 0 | 2 |
| d6 | single_joint_rate_MCC | 0.299980 | 0.587428 | 0.090250 | 3 | 33 | 0 | 14 |
| a4 | current | 0.522733 | 0.805784 | 0.465398 | 0 | 0 | 0 | 0 |
| a4 | single_selected_AP | 0.531759 | 0.810283 | 0.463866 | 74 | 74 | 17 | 18 |
| a4 | no_trunk_selected_AP | 0.528419 | 0.809880 | 0.460802 | 35 | 35 | 9 | 12 |
| a4 | multi_block_selected_AP | 0.527598 | 0.810128 | 0.440868 | 112 | 112 | 17 | 33 |
| a4 | multi_replace_selected_AP | 0.534418 | 0.809185 | 0.431656 | 153 | 153 | 26 | 48 |
| a4 | single_selected_MCC | 0.531162 | 0.809688 | 0.460802 | 61 | 61 | 13 | 16 |
| a4 | no_trunk_selected_MCC | 0.522733 | 0.805784 | 0.465398 | 0 | 0 | 0 | 0 |
| a4 | multi_block_selected_MCC | 0.522827 | 0.807571 | 0.453139 | 24 | 24 | 1 | 9 |
| a4 | multi_replace_selected_MCC | 0.531621 | 0.810365 | 0.462335 | 64 | 64 | 15 | 17 |
| a4 | single_joint_rate_MCC | 0.531162 | 0.809688 | 0.460802 | 61 | 61 | 13 | 16 |

All calls above rank labelled molecules only. Added/removed counts are diagnostic training calls, not blind-test changes. Raw probability Brier scores, fold-wise metrics, all fixed single placements and prior/rate curves are in the accompanying tables.

## Standalone representation screen

| target | scenario | AP | AUROC | MCC_target_prior |
| --- | --- | --- | --- | --- |
| d6 | standalone::core_separate_fe--trunk--lgbm_binary | 0.327746 | 0.616954 | 0.094847 |
| d6 | standalone::core_separate_hem_s--trunk--lgbm_binary | 0.313912 | 0.605478 | 0.094847 |
| d6 | standalone::core_mixed_hem--trunk--lgbm_binary | 0.308255 | 0.613799 | 0.084949 |
| d6 | standalone::core_separate_hem--trunk--lgbm_binary | 0.306625 | 0.587976 | 0.104965 |
| d6 | standalone::extended_separate_hem--trunk--lgbm_binary | 0.303645 | 0.590851 | 0.089871 |
| a4 | standalone::core_separate_fe--trunk_molecular--tabpfn_top500 | 0.500768 | 0.787917 | 0.388576 |
| a4 | standalone::ligand_only--trunk_molecular--tabpfn_top500 | 0.498256 | 0.787860 | 0.397821 |
| a4 | standalone::core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.496727 | 0.784551 | 0.382409 |
| a4 | standalone::core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.495956 | 0.790381 | 0.385493 |
| a4 | standalone::pxr_protein--trunk_molecular--tabpfn_top500 | 0.492517 | 0.782666 | 0.379324 |

These are global best-of-grid standalone readouts, shown only for representation diagnosis. They are not estimates of selected ensemble improvement.

## Paired conditional uncertainty

| target | scenario | reference | delta_AP | AP_lower | AP_upper | delta_MCC | MCC_lower | MCC_upper |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d6 | single_selected_AP | current | 0.009264 | 0.002484 | 0.016835 | -0.010231 | -0.030020 | 0.006037 |
| d6 | single_selected_AP | no_trunk_selected_AP | 0.012156 | 0.000654 | 0.023131 | -0.005087 | -0.028900 | 0.018059 |
| d6 | no_trunk_selected_AP | current | -0.002892 | -0.010904 | 0.006358 | -0.005144 | -0.032614 | 0.016870 |
| d6 | multi_block_selected_AP | current | 0.021949 | -0.008438 | 0.048718 | -0.005144 | -0.061440 | 0.042826 |
| d6 | multi_block_selected_AP | no_trunk_selected_AP | 0.024841 | -0.006417 | 0.052332 | 0.000000 | -0.054484 | 0.046261 |
| d6 | multi_replace_selected_AP | current | 0.003016 | -0.003580 | 0.010746 | -0.010231 | -0.031927 | 0.006585 |
| d6 | multi_replace_selected_AP | no_trunk_selected_AP | 0.005908 | -0.006990 | 0.018145 | -0.005087 | -0.032794 | 0.022490 |
| d6 | single_selected_MCC | current | 0.003332 | -0.005408 | 0.012342 | -0.010231 | -0.034054 | 0.010542 |
| d6 | single_selected_MCC | no_trunk_selected_MCC | 0.004970 | -0.005331 | 0.014998 | -0.010231 | -0.034809 | 0.010763 |
| d6 | no_trunk_selected_MCC | current | -0.001638 | -0.007647 | 0.004661 | 0.000000 | -0.013644 | 0.015086 |
| d6 | multi_block_selected_MCC | current | 0.005407 | -0.010899 | 0.025590 | -0.005144 | -0.032358 | 0.020531 |
| d6 | multi_block_selected_MCC | no_trunk_selected_MCC | 0.007045 | -0.010233 | 0.026494 | -0.005144 | -0.037476 | 0.024737 |
| d6 | multi_replace_selected_MCC | current | 0.001709 | -0.002265 | 0.006797 | -0.010231 | -0.025582 | 0.000000 |
| d6 | multi_replace_selected_MCC | no_trunk_selected_MCC | 0.003347 | -0.003857 | 0.010693 | -0.010231 | -0.031519 | 0.005878 |
| d6 | single_joint_rate_MCC | current | 0.003120 | -0.005955 | 0.011089 | -0.019859 | -0.055742 | 0.009628 |
| d6 | single_joint_rate_MCC | no_trunk_selected_MCC | 0.004758 | -0.006843 | 0.013668 | -0.019859 | -0.060622 | 0.011750 |
| a4 | single_selected_AP | current | 0.009027 | 0.000297 | 0.018518 | -0.001532 | -0.015694 | 0.012801 |
| a4 | single_selected_AP | no_trunk_selected_AP | 0.003340 | -0.003396 | 0.011720 | 0.003064 | -0.011686 | 0.017118 |
| a4 | no_trunk_selected_AP | current | 0.005687 | 0.002121 | 0.009288 | -0.004596 | -0.014983 | 0.006395 |
| a4 | multi_block_selected_AP | current | 0.004865 | -0.006505 | 0.018313 | -0.024530 | -0.045403 | -0.003445 |
| a4 | multi_block_selected_AP | no_trunk_selected_AP | -0.000822 | -0.011158 | 0.011234 | -0.019935 | -0.038967 | 0.000611 |
| a4 | multi_replace_selected_AP | current | 0.011685 | -0.002067 | 0.025060 | -0.033742 | -0.055335 | -0.010263 |
| a4 | multi_replace_selected_AP | no_trunk_selected_AP | 0.005998 | -0.005734 | 0.018089 | -0.029146 | -0.049799 | -0.007425 |
| a4 | single_selected_MCC | current | 0.008430 | 0.001155 | 0.015300 | -0.004596 | -0.017629 | 0.008801 |
| a4 | single_selected_MCC | no_trunk_selected_MCC | 0.008430 | 0.001484 | 0.014678 | -0.004596 | -0.017583 | 0.008456 |
| a4 | no_trunk_selected_MCC | current | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| a4 | multi_block_selected_MCC | current | 0.000094 | -0.009044 | 0.008452 | -0.012259 | -0.021583 | -0.003349 |
| a4 | multi_block_selected_MCC | no_trunk_selected_MCC | 0.000094 | -0.008518 | 0.008587 | -0.012259 | -0.021085 | -0.003530 |
| a4 | multi_replace_selected_MCC | current | 0.008889 | 0.001300 | 0.016237 | -0.003064 | -0.016845 | 0.010848 |
| a4 | multi_replace_selected_MCC | no_trunk_selected_MCC | 0.008889 | 0.001615 | 0.016586 | -0.003064 | -0.017048 | 0.010830 |
| a4 | single_joint_rate_MCC | current | 0.008430 | 0.001411 | 0.015389 | -0.004596 | -0.017621 | 0.009181 |
| a4 | single_joint_rate_MCC | no_trunk_selected_MCC | 0.008430 | 0.001393 | 0.015273 | -0.004596 | -0.017072 | 0.008593 |

2,000 paired UMAP-group resamples retain selected scores and calls. They do not rerank, refit selectors or account for searching the grid. Intervals are conditional descriptive uncertainty, not independent confirmation or multiplicity-adjusted evidence.

## Fold choices

| target | objective | family | fold | recipe | rate |
| --- | --- | --- | --- | --- | --- |
| d6 | AP | single | 0 | add::core_mixed_hem--trunk--lgbm_binary | 0.065333 |
| d6 | AP | no_trunk | 0 | original_add::pfn::external_mixed_pca256 | 0.065333 |
| d6 | AP | multi_block | 0 | {"kind": "blend", "count": 1, "weight": 0.5, "models": ["core_separate_fe--trunk--lgbm_binary"]} | 0.065333 |
| d6 | AP | multi_replace | 0 | [{"member": "pfn::external_mixed_pca256", "trunk": "core_separate_hem_s--trunk--lgbm_binary"}] | 0.065333 |
| d6 | AP | single | 1 | replace_pca_raw::core_mixed_hem--trunk--lgbm_binary | 0.065333 |
| d6 | AP | no_trunk | 1 | original_add::lgbm::external_mixed_plain | 0.065333 |
| d6 | AP | multi_block | 1 | {"kind": "blend", "count": 1, "weight": 0.5, "models": ["core_separate_fe--trunk--lgbm_binary"]} | 0.065333 |
| d6 | AP | multi_replace | 1 | [{"member": "pfn::external_mixed_pca256", "trunk": "core_mixed_hem--trunk--lgbm_binary"}] | 0.065333 |
| d6 | AP | single | 2 | add::core_mixed_hem--trunk--lgbm_binary | 0.065333 |
| d6 | AP | no_trunk | 2 | original_add::pfn::external_mixed_pca256 | 0.065333 |
| d6 | AP | multi_block | 2 | {"kind": "blend", "count": 5, "weight": 0.5, "models": ["core_separate_fe--trunk--lgbm_binary", "core_separate_fe--trunk_molecular--lgbm_binary", "ligand_only--trunk--tabpfn_top500", "pxr_core--trunk--lgbm_binary", "core_mixed_hem--trunk_molecular--lgbm_binary"]} | 0.065333 |
| d6 | AP | multi_replace | 2 | [{"member": "pfn::external_mixed_pca256", "trunk": "pxr_core--trunk--lgbm_binary"}] | 0.065333 |
| d6 | AP | single | 3 | add::core_mixed_hem--trunk--lgbm_binary | 0.065333 |
| d6 | AP | no_trunk | 3 | original_add::pfn::external_mixed_raw600 | 0.065333 |
| d6 | AP | multi_block | 3 | {"kind": "blend", "count": 1, "weight": 0.5, "models": ["core_separate_fe--trunk--lgbm_binary"]} | 0.065333 |
| d6 | AP | multi_replace | 3 | [{"member": "pfn::external_mixed_pca256", "trunk": "core_mixed_hem--trunk--lgbm_binary"}] | 0.065333 |
| d6 | AP | single | 4 | add::core_mixed_hem--trunk--lgbm_binary | 0.065333 |
| d6 | AP | no_trunk | 4 | original_add::pfn::external_mixed_pca256 | 0.065333 |
| d6 | AP | multi_block | 4 | {"kind": "blend", "count": 1, "weight": 0.5, "models": ["core_separate_fe--trunk--lgbm_binary"]} | 0.065333 |
| d6 | AP | multi_replace | 4 | [{"member": "pfn::external_mixed_pca256", "trunk": "core_separate_fe--trunk--lgbm_binary"}] | 0.065333 |
| d6 | MCC | single | 0 | add::core_distribution_hem--trunk--tabpfn_top500 | 0.065333 |
| d6 | MCC | no_trunk | 0 | original_add::graph::tdi_aux | 0.065333 |
| d6 | MCC | single_joint_rate | 0 | add::pxr_extended--trunk--lgbm_binary | 0.040000 |
| d6 | MCC | multi_block | 0 | {"kind": "blend", "count": 1, "weight": 0.25, "models": ["core_separate_hem--trunk--lgbm_binary"]} | 0.065333 |
| d6 | MCC | multi_replace | 0 | [{"member": "pfn::external_mixed_pca256", "trunk": "core_separate_hem_s--trunk--lgbm_binary"}] | 0.065333 |
| d6 | MCC | single | 1 | add::core_distribution_hem--trunk--tabpfn_top500 | 0.065333 |
| d6 | MCC | no_trunk | 1 | original_add::lgbm::morgan_plain | 0.065333 |
| d6 | MCC | single_joint_rate | 1 | replace_pca_raw::pxr_extended--trunk--lgbm_binary | 0.040000 |
| d6 | MCC | multi_block | 1 | {"kind": "current", "count": 0, "weight": 0.0, "models": []} | 0.065333 |
| d6 | MCC | multi_replace | 1 | [] | 0.065333 |
| d6 | MCC | single | 2 | add::core_distribution_hem--trunk--tabpfn_top500 | 0.065333 |
| d6 | MCC | no_trunk | 2 | original_add::graph::tdi_aux | 0.065333 |
| d6 | MCC | single_joint_rate | 2 | add::core_distribution_hem--trunk--tabpfn_top500 | 0.065333 |
| d6 | MCC | multi_block | 2 | {"kind": "blend", "count": 2, "weight": 0.25, "models": ["core_separate_hem--trunk--lgbm_binary", "extended_separate_hem--trunk_molecular--lgbm_binary"]} | 0.065333 |
| d6 | MCC | multi_replace | 2 | [] | 0.065333 |
| d6 | MCC | single | 3 | add::pxr_extended--trunk--lgbm_binary | 0.065333 |
| d6 | MCC | no_trunk | 3 | original_add::graph::tdi_aux | 0.065333 |
| d6 | MCC | single_joint_rate | 3 | add::extended_separate_hem--trunk--lgbm_binary | 0.040000 |
| d6 | MCC | multi_block | 3 | {"kind": "blend", "count": 1, "weight": 0.25, "models": ["core_separate_hem--trunk--lgbm_binary"]} | 0.065333 |
| d6 | MCC | multi_replace | 3 | [] | 0.065333 |
| d6 | MCC | single | 4 | add::pxr_core--trunk_molecular--tabpfn_top500 | 0.065333 |
| d6 | MCC | no_trunk | 4 | original_add::graph::tdi_aux | 0.065333 |
| d6 | MCC | single_joint_rate | 4 | add::pxr_extended--trunk--lgbm_binary | 0.040000 |
| d6 | MCC | multi_block | 4 | {"kind": "current", "count": 0, "weight": 0.0, "models": []} | 0.065333 |
| d6 | MCC | multi_replace | 4 | [] | 0.065333 |
| a4 | AP | single | 0 | replace_pxr_admetica::core_separate_hem_s--trunk--tabpfn_top500 | 0.332000 |
| a4 | AP | no_trunk | 0 | original_add::historical | 0.332000 |
| a4 | AP | multi_block | 0 | {"kind": "individual_add", "count": 5, "weight": 0.0, "models": ["core_separate_fe--trunk_molecular--tabpfn_top500", "ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500", "pxr_protein--trunk_molecular--lgbm_binary", "core_bidirectional_hem--trunk_molecular--tabpfn_top500"]} | 0.332000 |
| a4 | AP | multi_replace | 0 | [{"member": "pxr_admetica", "trunk": "core_separate_hem_s--trunk--tabpfn_top500"}, {"member": "jeremy_tabicl_pca256", "trunk": "core_bidirectional_hem--trunk_molecular--tabpfn_top500"}, {"member": "external_class_mordred", "trunk": "core_separate_hem_s--trunk_molecular--tabpfn_top500"}] | 0.332000 |
| a4 | AP | single | 1 | replace_pxr_admetica::core_bidirectional_hem--trunk--lgbm_binary | 0.332000 |
| a4 | AP | no_trunk | 1 | original_add::historical | 0.332000 |
| a4 | AP | multi_block | 1 | {"kind": "blend", "count": 1, "weight": 0.25, "models": ["core_separate_hem_s--trunk_molecular--tabpfn_top500"]} | 0.332000 |
| a4 | AP | multi_replace | 1 | [{"member": "pxr_admetica", "trunk": "core_bidirectional_hem--trunk--lgbm_binary"}, {"member": "jeremy_tabicl_pca256", "trunk": "core_separate_hem_s--trunk_molecular--tabpfn_top500"}, {"member": "external_class_mordred", "trunk": "core_bidirectional_hem--trunk_molecular--tabpfn_top500"}] | 0.332000 |
| a4 | AP | single | 2 | replace_pxr_admetica::core_separate_hem--trunk_molecular--tabpfn_top500 | 0.332000 |
| a4 | AP | no_trunk | 2 | original_add::historical | 0.332000 |
| a4 | AP | multi_block | 2 | {"kind": "individual_add", "count": 3, "weight": 0.0, "models": ["ligand_only--trunk_molecular--tabpfn_top500", "core_separate_fe--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500"]} | 0.332000 |
| a4 | AP | multi_replace | 2 | [{"member": "pxr_admetica", "trunk": "core_separate_hem--trunk_molecular--tabpfn_top500"}, {"member": "jeremy_tabicl_raw600", "trunk": "core_bidirectional_hem--trunk--lgbm_binary"}, {"member": "external_class_mordred", "trunk": "core_separate_hem_s--trunk_molecular--tabpfn_top500"}] | 0.332000 |
| a4 | AP | single | 3 | replace_pxr_admetica::core_separate_hem--trunk_molecular--tabpfn_top500 | 0.332000 |
| a4 | AP | no_trunk | 3 | original_add::historical | 0.332000 |
| a4 | AP | multi_block | 3 | {"kind": "blend", "count": 1, "weight": 0.5, "models": ["core_separate_fe--trunk_molecular--tabpfn_top500"]} | 0.332000 |
| a4 | AP | multi_replace | 3 | [{"member": "pxr_admetica", "trunk": "core_separate_hem--trunk_molecular--tabpfn_top500"}, {"member": "external_class_mordred", "trunk": "core_bidirectional_hem--trunk_molecular--lgbm_binary"}, {"member": "jeremy_tabicl_raw600", "trunk": "core_separate_hem_s--trunk_molecular--tabpfn_top500"}] | 0.332000 |
| a4 | AP | single | 4 | replace_pxr_admetica::core_separate_hem--trunk_molecular--tabpfn_top500 | 0.332000 |
| a4 | AP | no_trunk | 4 | original_add::historical | 0.332000 |
| a4 | AP | multi_block | 4 | {"kind": "blend", "count": 2, "weight": 0.5, "models": ["core_separate_fe--trunk_molecular--tabpfn_top500", "ligand_only--trunk_molecular--tabpfn_top500"]} | 0.332000 |
| a4 | AP | multi_replace | 4 | [{"member": "pxr_admetica", "trunk": "core_separate_hem--trunk_molecular--tabpfn_top500"}, {"member": "external_class_mordred", "trunk": "core_separate_fe--trunk_molecular--tabpfn_top500"}, {"member": "jeremy_tabicl_pca256", "trunk": "core_bidirectional_hem--trunk_molecular--tabpfn_top500"}] | 0.332000 |
| a4 | MCC | single | 0 | replace_jeremy_tabicl_raw600::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | no_trunk | 0 | current | 0.332000 |
| a4 | MCC | single_joint_rate | 0 | replace_jeremy_tabicl_raw600::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | multi_block | 0 | {"kind": "blend", "count": 1, "weight": 0.1, "models": ["ligand_only--trunk_molecular--tabpfn_top500"]} | 0.332000 |
| a4 | MCC | multi_replace | 0 | [{"member": "jeremy_tabicl_raw600", "trunk": "core_separate_fe--trunk--lgbm_binary"}] | 0.332000 |
| a4 | MCC | single | 1 | add::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | no_trunk | 1 | current | 0.332000 |
| a4 | MCC | single_joint_rate | 1 | add::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | multi_block | 1 | {"kind": "blend", "count": 5, "weight": 0.1, "models": ["ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500", "core_separate_hem_s--trunk--tabpfn_top500", "pxr_protein--trunk--tabpfn_top500", "core_separate_fe--trunk--tabpfn_top500"]} | 0.332000 |
| a4 | MCC | multi_replace | 1 | [{"member": "jeremy_tabicl_raw600", "trunk": "core_separate_fe--trunk--lgbm_binary"}] | 0.332000 |
| a4 | MCC | single | 2 | replace_jeremy_tabicl_raw600::core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.332000 |
| a4 | MCC | no_trunk | 2 | current | 0.332000 |
| a4 | MCC | single_joint_rate | 2 | replace_jeremy_tabicl_raw600::core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.332000 |
| a4 | MCC | multi_block | 2 | {"kind": "blend", "count": 1, "weight": 0.1, "models": ["pxr_extended--trunk_molecular--tabpfn_top500"]} | 0.332000 |
| a4 | MCC | multi_replace | 2 | [{"member": "jeremy_tabicl_raw600", "trunk": "core_mixed_hem--trunk_molecular--tabpfn_top500"}] | 0.332000 |
| a4 | MCC | single | 3 | replace_jeremy_tabicl_raw600::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | no_trunk | 3 | current | 0.332000 |
| a4 | MCC | single_joint_rate | 3 | replace_jeremy_tabicl_raw600::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | multi_block | 3 | {"kind": "blend", "count": 5, "weight": 0.25, "models": ["ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500", "core_bidirectional_hem--trunk_molecular--lgbm_binary", "extended_separate_hem--trunk--lgbm_binary", "pxr_extended--trunk_molecular--lgbm_binary"]} | 0.332000 |
| a4 | MCC | multi_replace | 3 | [{"member": "jeremy_tabicl_raw600", "trunk": "core_separate_fe--trunk--lgbm_binary"}] | 0.332000 |
| a4 | MCC | single | 4 | replace_jeremy_tabicl_raw600::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | no_trunk | 4 | current | 0.332000 |
| a4 | MCC | single_joint_rate | 4 | replace_jeremy_tabicl_raw600::core_separate_fe--trunk--lgbm_binary | 0.332000 |
| a4 | MCC | multi_block | 4 | {"kind": "current", "count": 0, "weight": 0.0, "models": []} | 0.332000 |
| a4 | MCC | multi_replace | 4 | [{"member": "jeremy_tabicl_raw600", "trunk": "core_separate_fe--trunk--lgbm_binary"}] | 0.332000 |

## Prediction and error complementarity

### d6: least-correlated readouts versus incumbent

| model | prediction_spearman | rank_residual_pearson | rank_squared_loss_pearson |
| --- | --- | --- | --- |
| pxr_protein--trunk--lgbm_binary | 0.401412 | 0.781727 | 0.364569 |
| pxr_core--trunk--lgbm_binary | 0.401855 | 0.777255 | 0.363140 |
| core_separate_hem--trunk--lgbm_binary | 0.414906 | 0.783585 | 0.378337 |
| core_separate_fe--trunk--lgbm_binary | 0.416743 | 0.779721 | 0.389259 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.417746 | 0.786039 | 0.396634 |
| extended_separate_hem--trunk--lgbm_binary | 0.418711 | 0.784548 | 0.375968 |
| core_separate_hem_s--trunk--lgbm_binary | 0.419217 | 0.782430 | 0.378448 |
| pxr_core--trunk--tabpfn_top500 | 0.423117 | 0.788751 | 0.382409 |

Trunk pair prediction Spearman range: 0.394898 to 0.884005; median 0.614532.

### a4: least-correlated readouts versus incumbent

| model | prediction_spearman | rank_residual_pearson | rank_squared_loss_pearson |
| --- | --- | --- | --- |
| core_mixed_hem--trunk--lgbm_binary | 0.651930 | 0.822807 | 0.590455 |
| pxr_core--trunk--lgbm_binary | 0.656763 | 0.824517 | 0.605217 |
| pxr_protein--trunk--lgbm_binary | 0.668128 | 0.829959 | 0.611559 |
| core_separate_hem_s--trunk--lgbm_binary | 0.674071 | 0.832972 | 0.612407 |
| extended_mixed_hem--trunk--lgbm_binary | 0.677629 | 0.833763 | 0.622406 |
| core_separate_hem--trunk--lgbm_binary | 0.678485 | 0.833919 | 0.615788 |
| core_separate_fe--trunk--lgbm_binary | 0.682997 | 0.835907 | 0.621991 |
| pxr_extended--trunk--lgbm_binary | 0.683609 | 0.838170 | 0.631314 |

Trunk pair prediction Spearman range: 0.714802 to 0.974607; median 0.834282.

Residuals here are within-fold percentile ranks minus binary labels. Their shared label term can produce high correlation, so low prediction correlation alone is not evidence of useful error diversity. Squared-rank-loss correlations and held-fold call swaps provide complementary diagnostics. These are not calibrated-probability residuals.

## Feature interpretation

| view | features |
| --- | --- |
| ligand_only | Ligand s mean; ligand-internal z mean/max |
| pxr_protein | Protein s mean + ligand s mean; whole-protein-to-ligand z mean/max |
| pxr_core | Protein s mean + ligand s mean; core-pocket-to-ligand z mean/max |
| pxr_extended | Protein s mean + ligand s mean; extended-pocket-to-ligand z mean/max |
| core_mixed_hem | Protein/ligand s means; core pocket and HEM pooled together for ligand z |
| extended_mixed_hem | Protein/ligand s means; extended pocket and HEM pooled together |
| core_separate_hem | Protein/ligand s means; core-pocket-to-ligand and HEM-to-ligand z separately |
| extended_separate_hem | Extended-pocket-to-ligand and HEM-to-ligand z separately |
| core_separate_fe | Core pocket, Fe and non-Fe HEM-to-ligand z separately |
| core_bidirectional_hem | Core pocket and HEM separate; both z directions to/from ligand |
| core_separate_hem_s | Core pocket/HEM separate z blocks, plus HEM s mean |
| core_distribution_hem | Core pocket/HEM separate z mean/std/q10/q90; ligand s mean/std |

`trunk_molecular` additionally uses frozen CheMeleon2048 + Mordred1613. `lgbm_binary` is the fixed LightGBM classifier; `tabpfn_top500` is TabPFN 2.6 with eight estimators and training-fold-only LightGBM gain selection up to 500 columns. Boltz z statistics are latent embeddings, not measured distances.

## Decision limits

A useful follow-up signal should improve the actual incumbent, survive held-fold selection and no-trunk composition controls, retain useful labelled positive additions at nearby rates, and show more than a global grid maximum. Fold instability, rate-grid boundary selection and wide paired intervals lower the case for storage-expensive full-pool extraction. AP increases alone do not establish MCC or leaderboard improvement.

The historical D6 AP-versus-control replay and A4 249-to-275 rate probe both failed to reproduce prospective leaderboard directions. Thus even an encouraging signal here is only grounds to consider completing the pool, not submission authorization.
