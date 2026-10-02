# CYP frozen Boltz trunk pooling screen

This is an exploratory fixed-grid comparison, not a submission decision. The Boltz trunk is frozen; all tensors include HEM. 'No HEM' variants omit HEM from readout pooling, not from the trunk input.

All models use the original five UMAP folds. Molecular controls use CheMeleon2048 + Mordred1613. Imputation and TabPFN top500 gain selection are fitted only on each training fold. Test predictions average five models.

Trunk features contain no predicted coordinates. Pair features are latent representations and are not interpreted as measured contacts.

Rankings and blend weights below are exploratory. No nested recipe selection or independent confirmation is claimed. Paired intervals are conditional on these fitted models, without multiplicity correction.

## cyp1a2 / regression

Run: `20260930T233551833106Z-e209b403`; launcher exit code: 1.

Export recovered without refitting in `20261001T004229979336Z-f08b913c`.

Completed recipes: 75. Recorded failures: 0.

| Recipe | ST-RAE | Spearman | Blend 10% | Blend 25% | Blend 50% |
|---|---:|---:|---:|---:|---:|
| core_bidirectional_hem--trunk--lgbm_l1 | 0.7757 | 0.5227 | 0.8155 | 0.7974 | 0.7763 |
| core_bidirectional_hem--trunk_molecular--lgbm_l1 | 0.7794 | 0.5292 | 0.8197 | 0.8064 | 0.7895 |
| core_separate_hem_s--trunk--lgbm_l1 | 0.7838 | 0.5217 | 0.8170 | 0.8007 | 0.7822 |
| core_separate_hem--trunk_molecular--lgbm_l1 | 0.7877 | 0.5223 | 0.8209 | 0.8094 | 0.7957 |
| extended_separate_hem--trunk_molecular--lgbm_l1 | 0.7877 | 0.5223 | 0.8209 | 0.8094 | 0.7957 |
| core_separate_hem_s--trunk_molecular--lgbm_l1 | 0.7892 | 0.5158 | 0.8205 | 0.8087 | 0.7949 |
| core_separate_fe--trunk_molecular--lgbm_l1 | 0.7893 | 0.5243 | 0.8208 | 0.8093 | 0.7949 |
| pxr_protein--trunk--lgbm_l1 | 0.7905 | 0.5062 | 0.8170 | 0.8012 | 0.7842 |
| ligand_only--trunk--tabpfn_top500 | 0.7908 | 0.5363 | 0.8356 | 0.8098 | 0.7821 |
| core_separate_hem--trunk--lgbm_l1 | 0.7917 | 0.5082 | 0.8177 | 0.8024 | 0.7852 |
| extended_separate_hem--trunk--lgbm_l1 | 0.7917 | 0.5082 | 0.8177 | 0.8024 | 0.7852 |
| core_distribution_hem--trunk--lgbm_l1 | 0.7919 | 0.5143 | 0.8182 | 0.8037 | 0.7879 |
| ligand_only--trunk_molecular--lgbm_l1 | 0.7925 | 0.5104 | 0.8230 | 0.8142 | 0.8022 |
| pxr_protein--trunk_molecular--lgbm_l1 | 0.7929 | 0.5184 | 0.8209 | 0.8096 | 0.7970 |
| core_distribution_hem--trunk_molecular--lgbm_l1 | 0.7941 | 0.5109 | 0.8212 | 0.8104 | 0.7979 |
| core_bidirectional_hem--trunk--tabpfn_top500 | 0.7946 | 0.5463 | 0.8340 | 0.8072 | 0.7801 |
| core_mixed_hem--trunk_molecular--lgbm_l1 | 0.7953 | 0.5053 | 0.8219 | 0.8119 | 0.7998 |
| extended_mixed_hem--trunk_molecular--lgbm_l1 | 0.7953 | 0.5053 | 0.8219 | 0.8119 | 0.7998 |
| core_separate_fe--trunk--lgbm_l1 | 0.7958 | 0.5113 | 0.8178 | 0.8029 | 0.7869 |
| pxr_core--trunk_molecular--lgbm_l1 | 0.7965 | 0.5011 | 0.8224 | 0.8129 | 0.8014 |
| pxr_extended--trunk_molecular--lgbm_l1 | 0.7965 | 0.5011 | 0.8224 | 0.8129 | 0.8014 |
| pxr_core--trunk--lgbm_l1 | 0.7981 | 0.4942 | 0.8189 | 0.8054 | 0.7905 |
| pxr_extended--trunk--lgbm_l1 | 0.7981 | 0.4942 | 0.8189 | 0.8054 | 0.7905 |
| core_mixed_hem--trunk--lgbm_l1 | 0.7987 | 0.5053 | 0.8183 | 0.8046 | 0.7904 |
| extended_mixed_hem--trunk--lgbm_l1 | 0.7987 | 0.5053 | 0.8183 | 0.8046 | 0.7904 |
| ligand_only--trunk--lgbm_l1 | 0.8014 | 0.4893 | 0.8197 | 0.8077 | 0.7953 |
| core_separate_hem_s--trunk--tabpfn_top500 | 0.8024 | 0.5272 | 0.8354 | 0.8106 | 0.7842 |
| core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.8058 | 0.5421 | 0.8404 | 0.8210 | 0.8012 |
| core_separate_fe--trunk--tabpfn_top500 | 0.8069 | 0.5246 | 0.8359 | 0.8117 | 0.7884 |
| pxr_core--trunk--tabpfn_top500 | 0.8085 | 0.5145 | 0.8365 | 0.8136 | 0.7908 |
| pxr_extended--trunk--tabpfn_top500 | 0.8085 | 0.5145 | 0.8365 | 0.8136 | 0.7908 |
| ligand_only--trunk_molecular--tabpfn_top500 | 0.8086 | 0.5308 | 0.8433 | 0.8267 | 0.8084 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.8093 | 0.5176 | 0.8371 | 0.8139 | 0.7911 |
| core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.8096 | 0.5462 | 0.8397 | 0.8197 | 0.8001 |
| core_separate_hem--trunk--tabpfn_top500 | 0.8120 | 0.5219 | 0.8365 | 0.8129 | 0.7891 |
| extended_separate_hem--trunk--tabpfn_top500 | 0.8120 | 0.5219 | 0.8365 | 0.8129 | 0.7891 |
| pxr_protein--trunk--tabpfn_top500 | 0.8124 | 0.5164 | 0.8369 | 0.8140 | 0.7924 |
| core_separate_fe--trunk_molecular--tabpfn_top500 | 0.8127 | 0.5332 | 0.8420 | 0.8244 | 0.8067 |
| core_mixed_hem--trunk--tabpfn_top500 | 0.8141 | 0.5256 | 0.8373 | 0.8144 | 0.7923 |
| extended_mixed_hem--trunk--tabpfn_top500 | 0.8141 | 0.5256 | 0.8373 | 0.8144 | 0.7923 |
| pxr_core--trunk_molecular--tabpfn_top500 | 0.8149 | 0.5263 | 0.8444 | 0.8298 | 0.8147 |
| pxr_extended--trunk_molecular--tabpfn_top500 | 0.8149 | 0.5263 | 0.8444 | 0.8298 | 0.8147 |
| core_separate_hem--trunk_molecular--lgbm_l2 | 0.8161 | 0.5146 | 0.8565 | 0.8385 | 0.8184 |
| extended_separate_hem--trunk_molecular--lgbm_l2 | 0.8161 | 0.5146 | 0.8565 | 0.8385 | 0.8184 |
| core_bidirectional_hem--trunk_molecular--lgbm_l2 | 0.8165 | 0.5165 | 0.8549 | 0.8351 | 0.8131 |
| core_bidirectional_hem--trunk--lgbm_l2 | 0.8173 | 0.5129 | 0.8508 | 0.8265 | 0.8019 |
| core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.8211 | 0.5264 | 0.8438 | 0.8291 | 0.8152 |
| extended_mixed_hem--trunk_molecular--tabpfn_top500 | 0.8211 | 0.5264 | 0.8438 | 0.8291 | 0.8152 |
| pxr_core--trunk_molecular--lgbm_l2 | 0.8216 | 0.5027 | 0.8573 | 0.8406 | 0.8224 |
| pxr_extended--trunk_molecular--lgbm_l2 | 0.8216 | 0.5027 | 0.8573 | 0.8406 | 0.8224 |
| core_separate_hem--trunk_molecular--tabpfn_top500 | 0.8218 | 0.5272 | 0.8423 | 0.8256 | 0.8104 |
| extended_separate_hem--trunk_molecular--tabpfn_top500 | 0.8218 | 0.5272 | 0.8423 | 0.8256 | 0.8104 |
| pxr_protein--trunk_molecular--lgbm_l2 | 0.8220 | 0.5068 | 0.8567 | 0.8387 | 0.8206 |
| pxr_protein--trunk_molecular--tabpfn_top500 | 0.8230 | 0.5253 | 0.8433 | 0.8281 | 0.8135 |
| pxr_protein--trunk--lgbm_l2 | 0.8245 | 0.4982 | 0.8506 | 0.8264 | 0.8028 |
| core_separate_fe--trunk_molecular--lgbm_l2 | 0.8259 | 0.5117 | 0.8568 | 0.8396 | 0.8215 |
| core_distribution_hem--trunk_molecular--tabpfn_top500 | 0.8276 | 0.5300 | 0.8426 | 0.8263 | 0.8119 |
| core_distribution_hem--trunk--lgbm_l2 | 0.8284 | 0.5004 | 0.8515 | 0.8284 | 0.8058 |
| core_distribution_hem--trunk_molecular--lgbm_l2 | 0.8288 | 0.4987 | 0.8572 | 0.8404 | 0.8237 |
| core_separate_hem_s--trunk_molecular--lgbm_l2 | 0.8299 | 0.5012 | 0.8572 | 0.8404 | 0.8234 |
| molecular_only--molecular--lgbm_l1 | 0.8301 | 0.4575 | 0.8301 | 0.8301 | 0.8301 |
| core_separate_fe--trunk--lgbm_l2 | 0.8316 | 0.4994 | 0.8518 | 0.8291 | 0.8079 |
| pxr_core--trunk--lgbm_l2 | 0.8348 | 0.4966 | 0.8520 | 0.8297 | 0.8085 |
| pxr_extended--trunk--lgbm_l2 | 0.8348 | 0.4966 | 0.8520 | 0.8297 | 0.8085 |
| ligand_only--trunk_molecular--lgbm_l2 | 0.8356 | 0.5042 | 0.8593 | 0.8457 | 0.8316 |
| core_separate_hem_s--trunk--lgbm_l2 | 0.8359 | 0.4951 | 0.8524 | 0.8306 | 0.8098 |
| core_separate_hem--trunk--lgbm_l2 | 0.8379 | 0.4974 | 0.8523 | 0.8307 | 0.8118 |
| extended_separate_hem--trunk--lgbm_l2 | 0.8379 | 0.4974 | 0.8523 | 0.8307 | 0.8118 |
| ligand_only--trunk--lgbm_l2 | 0.8389 | 0.4981 | 0.8538 | 0.8335 | 0.8164 |
| core_mixed_hem--trunk_molecular--lgbm_l2 | 0.8398 | 0.5024 | 0.8585 | 0.8433 | 0.8289 |
| extended_mixed_hem--trunk_molecular--lgbm_l2 | 0.8398 | 0.5024 | 0.8585 | 0.8433 | 0.8289 |
| core_mixed_hem--trunk--lgbm_l2 | 0.8474 | 0.4899 | 0.8508 | 0.8274 | 0.8073 |
| extended_mixed_hem--trunk--lgbm_l2 | 0.8474 | 0.4899 | 0.8508 | 0.8274 | 0.8073 |
| molecular_only--molecular--tabpfn_top500 | 0.8562 | 0.4776 | 0.8562 | 0.8562 | 0.8562 |
| molecular_only--molecular--lgbm_l2 | 0.8711 | 0.4583 | 0.8711 | 0.8711 | 0.8711 |

Paired group-bootstrap differences versus the same-readout molecular control are saved in `paired-uncertainty.parquet`. Negative ST-RAE differences favour the trunk recipe. The reference denominator is fixed at the complete OOF label mean; intervals describe grouped error differences at that denominator.


All fixed recipes are shown; the minimum/maximum is not an unbiased estimate of a selected model's future performance.

## cyp2c9 / regression

Run: `20261001T003655120573Z-71cfd9cb`; launcher exit code: 1.

Export recovered without refitting in `20261001T014152092601Z-b821c02a`.

Completed recipes: 75. Recorded failures: 0.

| Recipe | ST-RAE | Spearman | Blend 10% | Blend 25% | Blend 50% |
|---|---:|---:|---:|---:|---:|
| extended_mixed_hem--trunk_molecular--tabpfn_top500 | 0.6081 | 0.6670 | 0.6607 | 0.6371 | 0.6118 |
| extended_separate_hem--trunk_molecular--tabpfn_top500 | 0.6172 | 0.6656 | 0.6597 | 0.6365 | 0.6132 |
| pxr_core--trunk_molecular--tabpfn_top500 | 0.6219 | 0.6629 | 0.6619 | 0.6406 | 0.6186 |
| core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.6227 | 0.6668 | 0.6595 | 0.6359 | 0.6142 |
| pxr_extended--trunk_molecular--tabpfn_top500 | 0.6229 | 0.6599 | 0.6615 | 0.6398 | 0.6180 |
| core_separate_hem--trunk_molecular--tabpfn_top500 | 0.6283 | 0.6609 | 0.6619 | 0.6410 | 0.6213 |
| core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.6301 | 0.6527 | 0.6631 | 0.6440 | 0.6244 |
| core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.6311 | 0.6592 | 0.6600 | 0.6369 | 0.6173 |
| core_separate_fe--trunk_molecular--tabpfn_top500 | 0.6322 | 0.6563 | 0.6609 | 0.6394 | 0.6193 |
| core_distribution_hem--trunk_molecular--tabpfn_top500 | 0.6333 | 0.6613 | 0.6604 | 0.6383 | 0.6187 |
| extended_separate_hem--trunk_molecular--lgbm_l1 | 0.6337 | 0.6324 | 0.6692 | 0.6547 | 0.6383 |
| core_separate_hem_s--trunk_molecular--lgbm_l1 | 0.6356 | 0.6408 | 0.6696 | 0.6553 | 0.6401 |
| core_bidirectional_hem--trunk_molecular--lgbm_l1 | 0.6372 | 0.6341 | 0.6698 | 0.6559 | 0.6410 |
| core_separate_fe--trunk--tabpfn_top500 | 0.6380 | 0.6480 | 0.6530 | 0.6238 | 0.5991 |
| pxr_extended--trunk--tabpfn_top500 | 0.6380 | 0.6476 | 0.6547 | 0.6267 | 0.6033 |
| extended_mixed_hem--trunk--tabpfn_top500 | 0.6383 | 0.6442 | 0.6539 | 0.6252 | 0.6010 |
| extended_separate_hem--trunk--tabpfn_top500 | 0.6395 | 0.6496 | 0.6547 | 0.6269 | 0.6045 |
| pxr_protein--trunk_molecular--tabpfn_top500 | 0.6411 | 0.6566 | 0.6628 | 0.6438 | 0.6260 |
| core_separate_hem--trunk--tabpfn_top500 | 0.6415 | 0.6493 | 0.6538 | 0.6250 | 0.6019 |
| core_distribution_hem--trunk_molecular--lgbm_l1 | 0.6422 | 0.6309 | 0.6704 | 0.6575 | 0.6431 |
| core_mixed_hem--trunk_molecular--lgbm_l1 | 0.6426 | 0.6282 | 0.6715 | 0.6596 | 0.6467 |
| pxr_extended--trunk_molecular--lgbm_l1 | 0.6431 | 0.6336 | 0.6710 | 0.6588 | 0.6454 |
| pxr_protein--trunk_molecular--lgbm_l1 | 0.6434 | 0.6340 | 0.6707 | 0.6583 | 0.6447 |
| core_separate_fe--trunk_molecular--lgbm_l1 | 0.6436 | 0.6290 | 0.6706 | 0.6575 | 0.6433 |
| extended_mixed_hem--trunk_molecular--lgbm_l1 | 0.6441 | 0.6318 | 0.6715 | 0.6597 | 0.6468 |
| pxr_core--trunk--tabpfn_top500 | 0.6472 | 0.6440 | 0.6553 | 0.6280 | 0.6076 |
| pxr_protein--trunk--tabpfn_top500 | 0.6474 | 0.6458 | 0.6561 | 0.6298 | 0.6087 |
| core_bidirectional_hem--trunk--tabpfn_top500 | 0.6476 | 0.6458 | 0.6536 | 0.6254 | 0.6045 |
| core_separate_hem--trunk_molecular--lgbm_l1 | 0.6477 | 0.6310 | 0.6711 | 0.6593 | 0.6465 |
| pxr_protein--trunk_molecular--lgbm_l2 | 0.6489 | 0.6327 | 0.6896 | 0.6698 | 0.6492 |
| core_mixed_hem--trunk_molecular--lgbm_l2 | 0.6494 | 0.6415 | 0.6911 | 0.6732 | 0.6542 |
| core_bidirectional_hem--trunk_molecular--lgbm_l2 | 0.6507 | 0.6373 | 0.6893 | 0.6688 | 0.6477 |
| core_separate_hem--trunk_molecular--lgbm_l2 | 0.6520 | 0.6445 | 0.6895 | 0.6685 | 0.6476 |
| core_separate_hem_s--trunk_molecular--lgbm_l2 | 0.6537 | 0.6375 | 0.6895 | 0.6690 | 0.6488 |
| core_distribution_hem--trunk_molecular--lgbm_l2 | 0.6538 | 0.6356 | 0.6896 | 0.6696 | 0.6492 |
| core_separate_hem_s--trunk--tabpfn_top500 | 0.6543 | 0.6431 | 0.6548 | 0.6277 | 0.6074 |
| extended_separate_hem--trunk_molecular--lgbm_l2 | 0.6544 | 0.6394 | 0.6893 | 0.6687 | 0.6486 |
| core_mixed_hem--trunk--tabpfn_top500 | 0.6548 | 0.6343 | 0.6567 | 0.6309 | 0.6124 |
| ligand_only--trunk_molecular--tabpfn_top500 | 0.6551 | 0.6411 | 0.6688 | 0.6570 | 0.6462 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.6563 | 0.6411 | 0.6552 | 0.6279 | 0.6090 |
| pxr_core--trunk_molecular--lgbm_l1 | 0.6576 | 0.6228 | 0.6724 | 0.6624 | 0.6524 |
| pxr_extended--trunk_molecular--lgbm_l2 | 0.6596 | 0.6336 | 0.6911 | 0.6733 | 0.6542 |
| core_separate_fe--trunk_molecular--lgbm_l2 | 0.6628 | 0.6331 | 0.6907 | 0.6723 | 0.6549 |
| extended_mixed_hem--trunk_molecular--lgbm_l2 | 0.6635 | 0.6290 | 0.6921 | 0.6752 | 0.6581 |
| ligand_only--trunk_molecular--lgbm_l1 | 0.6664 | 0.6165 | 0.6759 | 0.6700 | 0.6641 |
| pxr_core--trunk_molecular--lgbm_l2 | 0.6692 | 0.6269 | 0.6926 | 0.6765 | 0.6608 |
| core_distribution_hem--trunk--lgbm_l1 | 0.6705 | 0.6233 | 0.6678 | 0.6541 | 0.6440 |
| core_bidirectional_hem--trunk--lgbm_l1 | 0.6733 | 0.6155 | 0.6675 | 0.6531 | 0.6434 |
| ligand_only--trunk_molecular--lgbm_l2 | 0.6758 | 0.6203 | 0.6965 | 0.6852 | 0.6736 |
| molecular_only--molecular--tabpfn_top500 | 0.6789 | 0.6151 | 0.6789 | 0.6789 | 0.6789 |
| pxr_protein--trunk--lgbm_l1 | 0.6792 | 0.6038 | 0.6676 | 0.6535 | 0.6453 |
| pxr_extended--trunk--lgbm_l1 | 0.6801 | 0.6014 | 0.6690 | 0.6562 | 0.6485 |
| molecular_only--molecular--lgbm_l1 | 0.6806 | 0.5934 | 0.6806 | 0.6806 | 0.6806 |
| extended_separate_hem--trunk--lgbm_l1 | 0.6815 | 0.6005 | 0.6684 | 0.6552 | 0.6471 |
| core_separate_fe--trunk--lgbm_l1 | 0.6830 | 0.6080 | 0.6694 | 0.6585 | 0.6526 |
| pxr_core--trunk--lgbm_l1 | 0.6851 | 0.6071 | 0.6695 | 0.6589 | 0.6543 |
| core_separate_hem--trunk--lgbm_l1 | 0.6853 | 0.6070 | 0.6690 | 0.6566 | 0.6508 |
| pxr_extended--trunk--lgbm_l2 | 0.6880 | 0.6042 | 0.6847 | 0.6616 | 0.6445 |
| ligand_only--trunk--tabpfn_top500 | 0.6888 | 0.6206 | 0.6620 | 0.6437 | 0.6341 |
| pxr_protein--trunk--lgbm_l2 | 0.6916 | 0.6059 | 0.6845 | 0.6616 | 0.6436 |
| core_separate_hem_s--trunk--lgbm_l1 | 0.6919 | 0.6023 | 0.6695 | 0.6587 | 0.6541 |
| core_separate_fe--trunk--lgbm_l2 | 0.6924 | 0.6011 | 0.6849 | 0.6626 | 0.6473 |
| core_distribution_hem--trunk--lgbm_l2 | 0.6933 | 0.6001 | 0.6845 | 0.6613 | 0.6453 |
| extended_mixed_hem--trunk--lgbm_l1 | 0.6936 | 0.5898 | 0.6698 | 0.6590 | 0.6550 |
| core_separate_hem--trunk--lgbm_l2 | 0.6941 | 0.6108 | 0.6846 | 0.6613 | 0.6452 |
| extended_separate_hem--trunk--lgbm_l2 | 0.6982 | 0.5940 | 0.6851 | 0.6638 | 0.6489 |
| core_bidirectional_hem--trunk--lgbm_l2 | 0.7003 | 0.6027 | 0.6860 | 0.6649 | 0.6511 |
| core_mixed_hem--trunk--lgbm_l1 | 0.7023 | 0.5858 | 0.6706 | 0.6609 | 0.6590 |
| pxr_core--trunk--lgbm_l2 | 0.7030 | 0.6070 | 0.6861 | 0.6655 | 0.6538 |
| extended_mixed_hem--trunk--lgbm_l2 | 0.7055 | 0.5870 | 0.6864 | 0.6661 | 0.6538 |
| molecular_only--molecular--lgbm_l2 | 0.7058 | 0.5870 | 0.7058 | 0.7058 | 0.7058 |
| core_mixed_hem--trunk--lgbm_l2 | 0.7109 | 0.5915 | 0.6865 | 0.6659 | 0.6532 |
| core_separate_hem_s--trunk--lgbm_l2 | 0.7135 | 0.5896 | 0.6872 | 0.6682 | 0.6589 |
| ligand_only--trunk--lgbm_l1 | 0.7203 | 0.5669 | 0.6739 | 0.6686 | 0.6731 |
| ligand_only--trunk--lgbm_l2 | 0.7317 | 0.5705 | 0.6917 | 0.6773 | 0.6733 |

Paired group-bootstrap differences versus the same-readout molecular control are saved in `paired-uncertainty.parquet`. Negative ST-RAE differences favour the trunk recipe. The reference denominator is fixed at the complete OOF label mean; intervals describe grouped error differences at that denominator.


All fixed recipes are shown; the minimum/maximum is not an unbiased estimate of a selected model's future performance.

## cyp2d6 / regression

Run: `20261001T013519598681Z-c88f83bd`; launcher exit code: 0.

Completed recipes: 75. Recorded failures: 0.

| Recipe | ST-RAE | Spearman | Blend 10% | Blend 25% | Blend 50% |
|---|---:|---:|---:|---:|---:|
| core_bidirectional_hem--trunk--lgbm_l1 | 0.8657 | 0.4302 | 0.9092 | 0.8937 | 0.8741 |
| core_separate_fe--trunk--lgbm_l1 | 0.8696 | 0.4199 | 0.9094 | 0.8945 | 0.8764 |
| core_distribution_hem--trunk_molecular--lgbm_l1 | 0.8724 | 0.4142 | 0.9123 | 0.9004 | 0.8846 |
| core_distribution_hem--trunk--lgbm_l1 | 0.8740 | 0.4124 | 0.9100 | 0.8957 | 0.8780 |
| core_separate_fe--trunk_molecular--lgbm_l1 | 0.8799 | 0.4117 | 0.9129 | 0.9022 | 0.8884 |
| extended_separate_hem--trunk--lgbm_l1 | 0.8805 | 0.4131 | 0.9106 | 0.8970 | 0.8810 |
| pxr_extended--trunk--lgbm_l1 | 0.8811 | 0.4090 | 0.9107 | 0.8973 | 0.8815 |
| extended_separate_hem--trunk_molecular--lgbm_l1 | 0.8812 | 0.4172 | 0.9135 | 0.9036 | 0.8910 |
| pxr_core--trunk_molecular--lgbm_l1 | 0.8819 | 0.4018 | 0.9134 | 0.9036 | 0.8908 |
| core_mixed_hem--trunk_molecular--lgbm_l1 | 0.8821 | 0.4051 | 0.9136 | 0.9036 | 0.8911 |
| core_separate_hem--trunk--lgbm_l1 | 0.8823 | 0.4073 | 0.9111 | 0.8982 | 0.8835 |
| core_bidirectional_hem--trunk_molecular--lgbm_l1 | 0.8835 | 0.4042 | 0.9132 | 0.9030 | 0.8897 |
| core_separate_hem--trunk_molecular--lgbm_l1 | 0.8840 | 0.4003 | 0.9136 | 0.9039 | 0.8916 |
| pxr_protein--trunk_molecular--lgbm_l1 | 0.8843 | 0.4100 | 0.9137 | 0.9044 | 0.8926 |
| core_mixed_hem--trunk--lgbm_l1 | 0.8846 | 0.4054 | 0.9109 | 0.8977 | 0.8833 |
| core_separate_hem_s--trunk_molecular--lgbm_l1 | 0.8851 | 0.3991 | 0.9139 | 0.9045 | 0.8924 |
| pxr_extended--trunk_molecular--lgbm_l1 | 0.8858 | 0.4106 | 0.9141 | 0.9051 | 0.8934 |
| extended_mixed_hem--trunk--lgbm_l1 | 0.8874 | 0.4049 | 0.9116 | 0.8994 | 0.8856 |
| core_separate_hem_s--trunk--lgbm_l1 | 0.8907 | 0.3969 | 0.9112 | 0.8988 | 0.8855 |
| extended_mixed_hem--trunk_molecular--lgbm_l1 | 0.8909 | 0.3994 | 0.9145 | 0.9062 | 0.8959 |
| pxr_core--trunk--lgbm_l1 | 0.8927 | 0.3883 | 0.9116 | 0.8999 | 0.8870 |
| pxr_protein--trunk--lgbm_l1 | 0.8955 | 0.3898 | 0.9122 | 0.9013 | 0.8896 |
| ligand_only--trunk_molecular--lgbm_l1 | 0.9014 | 0.3755 | 0.9172 | 0.9122 | 0.9054 |
| pxr_core--trunk--tabpfn_top500 | 0.9068 | 0.4136 | 0.9181 | 0.8982 | 0.8800 |
| ligand_only--trunk--lgbm_l1 | 0.9079 | 0.3432 | 0.9143 | 0.9061 | 0.8987 |
| core_separate_hem--trunk--tabpfn_top500 | 0.9088 | 0.4173 | 0.9182 | 0.8987 | 0.8813 |
| core_separate_hem_s--trunk--tabpfn_top500 | 0.9109 | 0.4225 | 0.9182 | 0.8986 | 0.8810 |
| core_separate_fe--trunk--tabpfn_top500 | 0.9128 | 0.4146 | 0.9190 | 0.9009 | 0.8857 |
| core_bidirectional_hem--trunk--tabpfn_top500 | 0.9142 | 0.4224 | 0.9191 | 0.9005 | 0.8845 |
| extended_separate_hem--trunk--tabpfn_top500 | 0.9176 | 0.4097 | 0.9184 | 0.8995 | 0.8835 |
| extended_mixed_hem--trunk--tabpfn_top500 | 0.9194 | 0.3990 | 0.9192 | 0.9015 | 0.8875 |
| molecular_only--molecular--lgbm_l1 | 0.9210 | 0.3471 | 0.9210 | 0.9210 | 0.9210 |
| core_separate_hem--trunk_molecular--tabpfn_top500 | 0.9219 | 0.4248 | 0.9245 | 0.9126 | 0.9024 |
| pxr_extended--trunk--tabpfn_top500 | 0.9222 | 0.4035 | 0.9192 | 0.9014 | 0.8874 |
| core_mixed_hem--trunk--tabpfn_top500 | 0.9225 | 0.4028 | 0.9198 | 0.9026 | 0.8884 |
| extended_separate_hem--trunk_molecular--tabpfn_top500 | 0.9268 | 0.4174 | 0.9257 | 0.9155 | 0.9069 |
| pxr_core--trunk_molecular--tabpfn_top500 | 0.9270 | 0.4173 | 0.9258 | 0.9156 | 0.9072 |
| core_separate_fe--trunk_molecular--tabpfn_top500 | 0.9273 | 0.4239 | 0.9254 | 0.9142 | 0.9063 |
| pxr_protein--trunk--tabpfn_top500 | 0.9278 | 0.3866 | 0.9204 | 0.9039 | 0.8922 |
| core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.9281 | 0.4137 | 0.9251 | 0.9137 | 0.9047 |
| core_distribution_hem--trunk_molecular--tabpfn_top500 | 0.9288 | 0.4159 | 0.9257 | 0.9154 | 0.9070 |
| ligand_only--trunk--tabpfn_top500 | 0.9297 | 0.3784 | 0.9216 | 0.9062 | 0.8941 |
| pxr_protein--trunk_molecular--tabpfn_top500 | 0.9298 | 0.4129 | 0.9269 | 0.9183 | 0.9117 |
| extended_mixed_hem--trunk_molecular--tabpfn_top500 | 0.9302 | 0.4093 | 0.9266 | 0.9177 | 0.9114 |
| core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.9315 | 0.4034 | 0.9266 | 0.9175 | 0.9100 |
| core_distribution_hem--trunk_molecular--lgbm_l2 | 0.9325 | 0.4081 | 0.9693 | 0.9512 | 0.9299 |
| pxr_extended--trunk_molecular--tabpfn_top500 | 0.9344 | 0.4143 | 0.9266 | 0.9175 | 0.9107 |
| molecular_only--molecular--tabpfn_top500 | 0.9349 | 0.4098 | 0.9349 | 0.9349 | 0.9349 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.9352 | 0.3936 | 0.9204 | 0.9044 | 0.8938 |
| ligand_only--trunk_molecular--tabpfn_top500 | 0.9397 | 0.4086 | 0.9307 | 0.9260 | 0.9234 |
| core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.9408 | 0.4132 | 0.9262 | 0.9168 | 0.9111 |
| extended_separate_hem--trunk_molecular--lgbm_l2 | 0.9420 | 0.4019 | 0.9708 | 0.9543 | 0.9366 |
| pxr_extended--trunk_molecular--lgbm_l2 | 0.9438 | 0.4016 | 0.9718 | 0.9568 | 0.9410 |
| core_bidirectional_hem--trunk_molecular--lgbm_l2 | 0.9458 | 0.3999 | 0.9697 | 0.9525 | 0.9363 |
| core_separate_hem--trunk_molecular--lgbm_l2 | 0.9498 | 0.3924 | 0.9707 | 0.9546 | 0.9392 |
| core_bidirectional_hem--trunk--lgbm_l2 | 0.9506 | 0.3965 | 0.9641 | 0.9405 | 0.9205 |
| pxr_core--trunk_molecular--lgbm_l2 | 0.9561 | 0.3837 | 0.9727 | 0.9596 | 0.9464 |
| core_separate_hem_s--trunk_molecular--lgbm_l2 | 0.9575 | 0.3844 | 0.9714 | 0.9563 | 0.9426 |
| core_separate_fe--trunk_molecular--lgbm_l2 | 0.9575 | 0.3884 | 0.9713 | 0.9561 | 0.9426 |
| core_distribution_hem--trunk--lgbm_l2 | 0.9577 | 0.3736 | 0.9642 | 0.9411 | 0.9216 |
| core_separate_fe--trunk--lgbm_l2 | 0.9618 | 0.3919 | 0.9647 | 0.9434 | 0.9262 |
| pxr_extended--trunk--lgbm_l2 | 0.9629 | 0.3823 | 0.9652 | 0.9432 | 0.9259 |
| extended_mixed_hem--trunk_molecular--lgbm_l2 | 0.9635 | 0.3822 | 0.9736 | 0.9610 | 0.9491 |
| core_separate_hem--trunk--lgbm_l2 | 0.9636 | 0.3902 | 0.9651 | 0.9439 | 0.9275 |
| pxr_protein--trunk_molecular--lgbm_l2 | 0.9648 | 0.3816 | 0.9739 | 0.9619 | 0.9509 |
| extended_separate_hem--trunk--lgbm_l2 | 0.9655 | 0.3926 | 0.9641 | 0.9415 | 0.9233 |
| core_mixed_hem--trunk_molecular--lgbm_l2 | 0.9665 | 0.3835 | 0.9737 | 0.9616 | 0.9506 |
| pxr_core--trunk--lgbm_l2 | 0.9668 | 0.3831 | 0.9665 | 0.9469 | 0.9313 |
| core_separate_hem_s--trunk--lgbm_l2 | 0.9689 | 0.3851 | 0.9646 | 0.9429 | 0.9261 |
| pxr_protein--trunk--lgbm_l2 | 0.9700 | 0.3775 | 0.9655 | 0.9441 | 0.9278 |
| ligand_only--trunk_molecular--lgbm_l2 | 0.9710 | 0.3613 | 0.9782 | 0.9716 | 0.9650 |
| ligand_only--trunk--lgbm_l2 | 0.9826 | 0.3333 | 0.9678 | 0.9504 | 0.9383 |
| molecular_only--molecular--lgbm_l2 | 0.9835 | 0.3457 | 0.9835 | 0.9835 | 0.9835 |
| extended_mixed_hem--trunk--lgbm_l2 | 0.9840 | 0.3631 | 0.9654 | 0.9445 | 0.9302 |
| core_mixed_hem--trunk--lgbm_l2 | 0.9945 | 0.3623 | 0.9675 | 0.9497 | 0.9408 |

Paired group-bootstrap differences versus the same-readout molecular control are saved in `paired-uncertainty.parquet`. Negative ST-RAE differences favour the trunk recipe. The reference denominator is fixed at the complete OOF label mean; intervals describe grouped error differences at that denominator.


All fixed recipes are shown; the minimum/maximum is not an unbiased estimate of a selected model's future performance.

## cyp2d6 / classification

Run: `20261001T023655009132Z-91096754`; launcher exit code: 0.

Completed recipes: 50. Recorded failures: 0.

TDI is a labelled-only diagnostic within the target-extracted pool. The complete 6,146-molecule operating pool is absent, so AP/AUROC and MCC at 0.5 do not support call-rate selection or deployment.

| Recipe | AP | AUROC | MCC at 0.5 |
|---|---:|---:|---:|
| core_separate_fe--trunk--lgbm_binary | 0.3245 | 0.6139 | 0.1137 |
| core_separate_hem_s--trunk--lgbm_binary | 0.3137 | 0.6016 | 0.1274 |
| core_separate_hem--trunk--lgbm_binary | 0.3108 | 0.5852 | 0.1322 |
| extended_separate_hem--trunk--lgbm_binary | 0.3102 | 0.5918 | 0.1028 |
| core_mixed_hem--trunk--lgbm_binary | 0.3064 | 0.6120 | 0.0872 |
| core_separate_fe--trunk_molecular--lgbm_binary | 0.3039 | 0.5952 | 0.0501 |
| pxr_extended--trunk--lgbm_binary | 0.3017 | 0.5950 | 0.1119 |
| pxr_core--trunk--lgbm_binary | 0.3017 | 0.5956 | 0.0717 |
| core_mixed_hem--trunk_molecular--lgbm_binary | 0.2996 | 0.5840 | 0.0732 |
| pxr_core--trunk_molecular--lgbm_binary | 0.2938 | 0.5792 | 0.0400 |
| core_distribution_hem--trunk_molecular--tabpfn_top500 | 0.2880 | 0.5848 | 0.0875 |
| core_separate_hem_s--trunk_molecular--lgbm_binary | 0.2864 | 0.5696 | 0.0382 |
| extended_separate_hem--trunk_molecular--lgbm_binary | 0.2845 | 0.5837 | -0.0005 |
| core_separate_fe--trunk--tabpfn_top500 | 0.2840 | 0.5771 | 0.0377 |
| pxr_protein--trunk_molecular--lgbm_binary | 0.2828 | 0.5665 | 0.0508 |
| extended_mixed_hem--trunk--lgbm_binary | 0.2821 | 0.5949 | 0.0489 |
| core_bidirectional_hem--trunk--lgbm_binary | 0.2818 | 0.5806 | 0.0519 |
| core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.2815 | 0.5814 | 0.0141 |
| pxr_core--trunk--tabpfn_top500 | 0.2815 | 0.5760 | 0.0042 |
| core_mixed_hem--trunk--tabpfn_top500 | 0.2810 | 0.5944 | -0.0333 |
| core_separate_hem--trunk_molecular--lgbm_binary | 0.2799 | 0.5656 | 0.0557 |
| pxr_extended--trunk_molecular--lgbm_binary | 0.2799 | 0.5711 | 0.0437 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.2797 | 0.5764 | 0.0752 |
| core_distribution_hem--trunk--lgbm_binary | 0.2783 | 0.5740 | 0.1008 |
| pxr_extended--trunk--tabpfn_top500 | 0.2774 | 0.5801 | 0.0042 |
| extended_mixed_hem--trunk--tabpfn_top500 | 0.2770 | 0.5860 | 0.0042 |
| extended_mixed_hem--trunk_molecular--lgbm_binary | 0.2768 | 0.5698 | 0.0448 |
| pxr_core--trunk_molecular--tabpfn_top500 | 0.2763 | 0.5752 | 0.0543 |
| pxr_protein--trunk--lgbm_binary | 0.2761 | 0.5691 | 0.0472 |
| pxr_protein--trunk--tabpfn_top500 | 0.2756 | 0.5845 | 0.0127 |
| core_bidirectional_hem--trunk_molecular--lgbm_binary | 0.2747 | 0.5720 | 0.0497 |
| core_separate_hem--trunk_molecular--tabpfn_top500 | 0.2746 | 0.5670 | 0.0328 |
| extended_separate_hem--trunk--tabpfn_top500 | 0.2741 | 0.5849 | 0.0115 |
| core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.2731 | 0.5797 | 0.0290 |
| core_distribution_hem--trunk_molecular--lgbm_binary | 0.2727 | 0.5569 | 0.0815 |
| ligand_only--trunk--tabpfn_top500 | 0.2726 | 0.5643 | -0.0023 |
| core_separate_hem--trunk--tabpfn_top500 | 0.2713 | 0.5776 | 0.0011 |
| extended_mixed_hem--trunk_molecular--tabpfn_top500 | 0.2710 | 0.5720 | 0.0522 |
| ligand_only--trunk--lgbm_binary | 0.2709 | 0.5750 | 0.0135 |
| pxr_extended--trunk_molecular--tabpfn_top500 | 0.2705 | 0.5656 | 0.0454 |
| ligand_only--trunk_molecular--lgbm_binary | 0.2699 | 0.5687 | 0.0118 |
| ligand_only--trunk_molecular--tabpfn_top500 | 0.2689 | 0.5667 | 0.0454 |
| core_separate_fe--trunk_molecular--tabpfn_top500 | 0.2670 | 0.5708 | 0.0472 |
| pxr_protein--trunk_molecular--tabpfn_top500 | 0.2662 | 0.5668 | 0.0581 |
| core_separate_hem_s--trunk--tabpfn_top500 | 0.2660 | 0.5611 | 0.0060 |
| core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.2648 | 0.5625 | 0.0263 |
| extended_separate_hem--trunk_molecular--tabpfn_top500 | 0.2628 | 0.5688 | 0.0205 |
| molecular_only--molecular--lgbm_binary | 0.2601 | 0.5476 | 0.0565 |
| core_bidirectional_hem--trunk--tabpfn_top500 | 0.2582 | 0.5666 | 0.0242 |
| molecular_only--molecular--tabpfn_top500 | 0.2433 | 0.5319 | -0.0072 |

All fixed recipes are shown; the minimum/maximum is not an unbiased estimate of a selected model's future performance.

## cyp3a4 / regression

Run: `20261001T031820053553Z-73768adb`; launcher exit code: 0.

Completed recipes: 75. Recorded failures: 0.

| Recipe | ST-RAE | Spearman | Blend 10% | Blend 25% | Blend 50% |
|---|---:|---:|---:|---:|---:|
| pxr_core--trunk_molecular--tabpfn_top500 | 0.5211 | 0.7757 | 0.5549 | 0.5407 | 0.5246 |
| ligand_only--trunk_molecular--tabpfn_top500 | 0.5221 | 0.7740 | 0.5558 | 0.5425 | 0.5275 |
| pxr_protein--trunk_molecular--tabpfn_top500 | 0.5230 | 0.7742 | 0.5555 | 0.5419 | 0.5265 |
| pxr_extended--trunk_molecular--tabpfn_top500 | 0.5231 | 0.7759 | 0.5551 | 0.5408 | 0.5248 |
| core_separate_fe--trunk_molecular--tabpfn_top500 | 0.5265 | 0.7730 | 0.5550 | 0.5409 | 0.5257 |
| core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.5290 | 0.7738 | 0.5540 | 0.5392 | 0.5243 |
| core_distribution_hem--trunk_molecular--tabpfn_top500 | 0.5292 | 0.7726 | 0.5543 | 0.5398 | 0.5247 |
| extended_mixed_hem--trunk_molecular--tabpfn_top500 | 0.5292 | 0.7716 | 0.5557 | 0.5424 | 0.5286 |
| core_separate_hem--trunk_molecular--tabpfn_top500 | 0.5300 | 0.7729 | 0.5557 | 0.5425 | 0.5287 |
| extended_separate_hem--trunk_molecular--tabpfn_top500 | 0.5319 | 0.7733 | 0.5547 | 0.5404 | 0.5260 |
| core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.5344 | 0.7689 | 0.5567 | 0.5446 | 0.5319 |
| core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.5345 | 0.7706 | 0.5557 | 0.5425 | 0.5288 |
| pxr_extended--trunk_molecular--lgbm_l1 | 0.5371 | 0.7501 | 0.5703 | 0.5591 | 0.5455 |
| pxr_extended--trunk_molecular--lgbm_l2 | 0.5417 | 0.7596 | 0.5804 | 0.5658 | 0.5484 |
| extended_separate_hem--trunk_molecular--lgbm_l2 | 0.5423 | 0.7575 | 0.5804 | 0.5659 | 0.5486 |
| core_separate_fe--trunk_molecular--lgbm_l1 | 0.5451 | 0.7451 | 0.5716 | 0.5620 | 0.5507 |
| extended_separate_hem--trunk_molecular--lgbm_l1 | 0.5452 | 0.7489 | 0.5710 | 0.5610 | 0.5496 |
| core_separate_hem--trunk_molecular--lgbm_l1 | 0.5465 | 0.7442 | 0.5718 | 0.5626 | 0.5520 |
| core_distribution_hem--trunk_molecular--lgbm_l1 | 0.5477 | 0.7438 | 0.5715 | 0.5620 | 0.5510 |
| core_bidirectional_hem--trunk_molecular--lgbm_l2 | 0.5490 | 0.7539 | 0.5813 | 0.5684 | 0.5532 |
| extended_mixed_hem--trunk_molecular--lgbm_l1 | 0.5498 | 0.7421 | 0.5722 | 0.5636 | 0.5537 |
| core_separate_hem_s--trunk_molecular--lgbm_l2 | 0.5508 | 0.7536 | 0.5816 | 0.5690 | 0.5541 |
| core_mixed_hem--trunk_molecular--lgbm_l1 | 0.5511 | 0.7444 | 0.5727 | 0.5646 | 0.5554 |
| core_separate_hem--trunk_molecular--lgbm_l2 | 0.5513 | 0.7538 | 0.5819 | 0.5697 | 0.5554 |
| pxr_protein--trunk_molecular--lgbm_l1 | 0.5516 | 0.7444 | 0.5728 | 0.5649 | 0.5558 |
| core_separate_hem_s--trunk_molecular--lgbm_l1 | 0.5523 | 0.7427 | 0.5725 | 0.5643 | 0.5550 |
| pxr_core--trunk_molecular--lgbm_l2 | 0.5525 | 0.7544 | 0.5825 | 0.5709 | 0.5572 |
| extended_separate_hem--trunk--tabpfn_top500 | 0.5531 | 0.7564 | 0.5499 | 0.5315 | 0.5174 |
| ligand_only--trunk_molecular--lgbm_l2 | 0.5544 | 0.7472 | 0.5840 | 0.5740 | 0.5619 |
| pxr_core--trunk_molecular--lgbm_l1 | 0.5553 | 0.7424 | 0.5733 | 0.5659 | 0.5577 |
| core_mixed_hem--trunk_molecular--lgbm_l2 | 0.5554 | 0.7478 | 0.5833 | 0.5725 | 0.5601 |
| core_separate_fe--trunk_molecular--lgbm_l2 | 0.5570 | 0.7504 | 0.5824 | 0.5708 | 0.5578 |
| pxr_protein--trunk--tabpfn_top500 | 0.5579 | 0.7524 | 0.5497 | 0.5309 | 0.5164 |
| core_distribution_hem--trunk_molecular--lgbm_l2 | 0.5582 | 0.7480 | 0.5826 | 0.5713 | 0.5589 |
| core_bidirectional_hem--trunk_molecular--lgbm_l1 | 0.5588 | 0.7421 | 0.5727 | 0.5651 | 0.5577 |
| ligand_only--trunk_molecular--lgbm_l1 | 0.5592 | 0.7343 | 0.5745 | 0.5686 | 0.5621 |
| pxr_core--trunk--tabpfn_top500 | 0.5594 | 0.7530 | 0.5504 | 0.5322 | 0.5189 |
| core_bidirectional_hem--trunk--tabpfn_top500 | 0.5595 | 0.7532 | 0.5502 | 0.5324 | 0.5186 |
| pxr_protein--trunk_molecular--lgbm_l2 | 0.5596 | 0.7471 | 0.5837 | 0.5736 | 0.5621 |
| core_separate_hem_s--trunk--tabpfn_top500 | 0.5602 | 0.7519 | 0.5506 | 0.5330 | 0.5206 |
| pxr_extended--trunk--tabpfn_top500 | 0.5604 | 0.7532 | 0.5501 | 0.5319 | 0.5187 |
| core_separate_hem--trunk--tabpfn_top500 | 0.5606 | 0.7521 | 0.5506 | 0.5333 | 0.5212 |
| core_separate_fe--trunk--tabpfn_top500 | 0.5606 | 0.7533 | 0.5506 | 0.5335 | 0.5211 |
| extended_mixed_hem--trunk_molecular--lgbm_l2 | 0.5607 | 0.7448 | 0.5836 | 0.5733 | 0.5613 |
| molecular_only--molecular--tabpfn_top500 | 0.5665 | 0.7438 | 0.5665 | 0.5665 | 0.5665 |
| extended_mixed_hem--trunk--tabpfn_top500 | 0.5698 | 0.7416 | 0.5510 | 0.5345 | 0.5251 |
| ligand_only--trunk--tabpfn_top500 | 0.5699 | 0.7399 | 0.5518 | 0.5367 | 0.5277 |
| molecular_only--molecular--lgbm_l1 | 0.5792 | 0.7164 | 0.5792 | 0.5792 | 0.5792 |
| core_mixed_hem--trunk--tabpfn_top500 | 0.5799 | 0.7377 | 0.5515 | 0.5364 | 0.5293 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.5813 | 0.7465 | 0.5520 | 0.5372 | 0.5303 |
| core_distribution_hem--trunk--lgbm_l1 | 0.5823 | 0.7195 | 0.5671 | 0.5539 | 0.5456 |
| extended_separate_hem--trunk--lgbm_l2 | 0.5837 | 0.7254 | 0.5747 | 0.5554 | 0.5422 |
| core_bidirectional_hem--trunk--lgbm_l2 | 0.5848 | 0.7262 | 0.5756 | 0.5576 | 0.5446 |
| extended_separate_hem--trunk--lgbm_l1 | 0.5850 | 0.7210 | 0.5669 | 0.5535 | 0.5452 |
| core_bidirectional_hem--trunk--lgbm_l1 | 0.5858 | 0.7199 | 0.5680 | 0.5559 | 0.5483 |
| core_separate_hem--trunk--lgbm_l2 | 0.5868 | 0.7262 | 0.5757 | 0.5579 | 0.5457 |
| pxr_extended--trunk--lgbm_l2 | 0.5874 | 0.7258 | 0.5748 | 0.5555 | 0.5426 |
| pxr_core--trunk--lgbm_l2 | 0.5875 | 0.7237 | 0.5757 | 0.5579 | 0.5456 |
| pxr_extended--trunk--lgbm_l1 | 0.5878 | 0.7162 | 0.5674 | 0.5550 | 0.5469 |
| core_separate_fe--trunk--lgbm_l2 | 0.5880 | 0.7221 | 0.5756 | 0.5577 | 0.5451 |
| core_separate_hem--trunk--lgbm_l1 | 0.5885 | 0.7165 | 0.5677 | 0.5554 | 0.5482 |
| core_separate_hem_s--trunk--lgbm_l2 | 0.5907 | 0.7219 | 0.5756 | 0.5577 | 0.5472 |
| pxr_protein--trunk--lgbm_l1 | 0.5907 | 0.7122 | 0.5691 | 0.5583 | 0.5516 |
| core_separate_fe--trunk--lgbm_l1 | 0.5911 | 0.7167 | 0.5679 | 0.5561 | 0.5485 |
| molecular_only--molecular--lgbm_l2 | 0.5916 | 0.7218 | 0.5916 | 0.5916 | 0.5916 |
| core_separate_hem_s--trunk--lgbm_l1 | 0.5923 | 0.7155 | 0.5681 | 0.5564 | 0.5506 |
| core_distribution_hem--trunk--lgbm_l2 | 0.5928 | 0.7193 | 0.5768 | 0.5602 | 0.5497 |
| pxr_core--trunk--lgbm_l1 | 0.5945 | 0.7115 | 0.5676 | 0.5556 | 0.5498 |
| extended_mixed_hem--trunk--lgbm_l1 | 0.6045 | 0.6952 | 0.5696 | 0.5598 | 0.5558 |
| ligand_only--trunk--lgbm_l2 | 0.6068 | 0.7039 | 0.5781 | 0.5635 | 0.5563 |
| pxr_protein--trunk--lgbm_l2 | 0.6069 | 0.7104 | 0.5776 | 0.5624 | 0.5551 |
| ligand_only--trunk--lgbm_l1 | 0.6099 | 0.6864 | 0.5702 | 0.5620 | 0.5612 |
| extended_mixed_hem--trunk--lgbm_l2 | 0.6158 | 0.7024 | 0.5785 | 0.5652 | 0.5597 |
| core_mixed_hem--trunk--lgbm_l2 | 0.6242 | 0.6922 | 0.5790 | 0.5662 | 0.5631 |
| core_mixed_hem--trunk--lgbm_l1 | 0.6253 | 0.6802 | 0.5709 | 0.5641 | 0.5663 |

Paired group-bootstrap differences versus the same-readout molecular control are saved in `paired-uncertainty.parquet`. Negative ST-RAE differences favour the trunk recipe. The reference denominator is fixed at the complete OOF label mean; intervals describe grouped error differences at that denominator.


All fixed recipes are shown; the minimum/maximum is not an unbiased estimate of a selected model's future performance.

## cyp3a4 / classification

Run: `20261001T043809889682Z-cd81b014`; launcher exit code: 0.

Completed recipes: 50. Recorded failures: 0.

TDI is a labelled-only diagnostic within the target-extracted pool. The complete 6,146-molecule operating pool is absent, so AP/AUROC and MCC at 0.5 do not support call-rate selection or deployment.

| Recipe | AP | AUROC | MCC at 0.5 |
|---|---:|---:|---:|
| core_separate_fe--trunk_molecular--tabpfn_top500 | 0.4976 | 0.7887 | 0.3160 |
| core_separate_hem_s--trunk_molecular--tabpfn_top500 | 0.4949 | 0.7859 | 0.2953 |
| core_bidirectional_hem--trunk_molecular--tabpfn_top500 | 0.4935 | 0.7924 | 0.3114 |
| pxr_protein--trunk_molecular--tabpfn_top500 | 0.4916 | 0.7842 | 0.2850 |
| ligand_only--trunk_molecular--tabpfn_top500 | 0.4912 | 0.7867 | 0.2808 |
| core_mixed_hem--trunk_molecular--tabpfn_top500 | 0.4903 | 0.7853 | 0.2974 |
| core_separate_hem--trunk_molecular--tabpfn_top500 | 0.4875 | 0.7852 | 0.2900 |
| pxr_extended--trunk_molecular--tabpfn_top500 | 0.4871 | 0.7836 | 0.2715 |
| extended_mixed_hem--trunk_molecular--tabpfn_top500 | 0.4861 | 0.7822 | 0.2858 |
| pxr_core--trunk_molecular--tabpfn_top500 | 0.4845 | 0.7852 | 0.2808 |
| extended_separate_hem--trunk_molecular--tabpfn_top500 | 0.4821 | 0.7799 | 0.2921 |
| molecular_only--molecular--tabpfn_top500 | 0.4816 | 0.7797 | 0.2680 |
| extended_mixed_hem--trunk_molecular--lgbm_binary | 0.4810 | 0.7692 | 0.2393 |
| core_distribution_hem--trunk_molecular--tabpfn_top500 | 0.4772 | 0.7792 | 0.3103 |
| extended_separate_hem--trunk_molecular--lgbm_binary | 0.4764 | 0.7672 | 0.2402 |
| pxr_protein--trunk_molecular--lgbm_binary | 0.4753 | 0.7637 | 0.2429 |
| pxr_extended--trunk_molecular--lgbm_binary | 0.4739 | 0.7690 | 0.2472 |
| core_bidirectional_hem--trunk_molecular--lgbm_binary | 0.4735 | 0.7714 | 0.2453 |
| core_separate_hem--trunk_molecular--lgbm_binary | 0.4712 | 0.7649 | 0.2683 |
| core_separate_hem_s--trunk_molecular--lgbm_binary | 0.4700 | 0.7650 | 0.2548 |
| core_mixed_hem--trunk_molecular--lgbm_binary | 0.4695 | 0.7620 | 0.2322 |
| ligand_only--trunk_molecular--lgbm_binary | 0.4672 | 0.7618 | 0.2503 |
| core_distribution_hem--trunk_molecular--lgbm_binary | 0.4656 | 0.7630 | 0.2412 |
| core_separate_fe--trunk_molecular--lgbm_binary | 0.4653 | 0.7625 | 0.2541 |
| pxr_core--trunk_molecular--lgbm_binary | 0.4607 | 0.7585 | 0.2346 |
| core_bidirectional_hem--trunk--tabpfn_top500 | 0.4602 | 0.7716 | 0.2341 |
| extended_mixed_hem--trunk--tabpfn_top500 | 0.4572 | 0.7646 | 0.2356 |
| pxr_extended--trunk--tabpfn_top500 | 0.4555 | 0.7616 | 0.2541 |
| extended_separate_hem--trunk--tabpfn_top500 | 0.4550 | 0.7653 | 0.2521 |
| molecular_only--molecular--lgbm_binary | 0.4522 | 0.7535 | 0.2263 |
| core_distribution_hem--trunk--tabpfn_top500 | 0.4515 | 0.7569 | 0.2640 |
| ligand_only--trunk--tabpfn_top500 | 0.4512 | 0.7539 | 0.2531 |
| core_separate_fe--trunk--tabpfn_top500 | 0.4449 | 0.7614 | 0.2578 |
| core_separate_hem--trunk--tabpfn_top500 | 0.4432 | 0.7626 | 0.2124 |
| core_separate_hem_s--trunk--tabpfn_top500 | 0.4410 | 0.7644 | 0.2074 |
| core_bidirectional_hem--trunk--lgbm_binary | 0.4393 | 0.7565 | 0.2134 |
| pxr_protein--trunk--tabpfn_top500 | 0.4356 | 0.7570 | 0.2141 |
| extended_mixed_hem--trunk--lgbm_binary | 0.4321 | 0.7436 | 0.2034 |
| core_mixed_hem--trunk--tabpfn_top500 | 0.4317 | 0.7542 | 0.1950 |
| pxr_core--trunk--tabpfn_top500 | 0.4295 | 0.7507 | 0.2110 |
| ligand_only--trunk--lgbm_binary | 0.4291 | 0.7380 | 0.2090 |
| pxr_extended--trunk--lgbm_binary | 0.4282 | 0.7382 | 0.2275 |
| pxr_protein--trunk--lgbm_binary | 0.4280 | 0.7394 | 0.2350 |
| core_separate_fe--trunk--lgbm_binary | 0.4261 | 0.7470 | 0.2181 |
| extended_separate_hem--trunk--lgbm_binary | 0.4260 | 0.7464 | 0.2097 |
| core_separate_hem_s--trunk--lgbm_binary | 0.4186 | 0.7397 | 0.1924 |
| core_separate_hem--trunk--lgbm_binary | 0.4135 | 0.7454 | 0.1947 |
| core_mixed_hem--trunk--lgbm_binary | 0.4126 | 0.7332 | 0.1639 |
| core_distribution_hem--trunk--lgbm_binary | 0.4119 | 0.7427 | 0.1929 |
| pxr_core--trunk--lgbm_binary | 0.4098 | 0.7367 | 0.1969 |

All fixed recipes are shown; the minimum/maximum is not an unbiased estimate of a selected model's future performance.

## Failures and scope

[]

Original large tensor and pooled feature payloads stay local with recorded input hashes. Completed runs and negative results are retained. No candidate or submission is uploaded.
