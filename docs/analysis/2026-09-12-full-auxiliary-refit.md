# Full-data auxiliary refit for the five-component submission

Run `shared/20260912T074339967571Z-ff194fcf` refits the adopted Rustuna-selected
CheMeleon log2fc recipe on all 4,376 compounds with measured primary-screen
responses (17,504 labels, four CYP endpoints). No validation subset is retained
in this deployment fit. Training lasts 13 epochs, the rounded median of the ten
existing auxiliary pair models' selected epoch counts. The original 80-epoch
learning-rate decay horizon is retained so that the selected epoch has the same
learning-rate position. Target centering/scaling uses all measured training labels.

`outputs/full-aux.parquet` contains 6,896 unique compound IDs and four feature
columns `aux_0` through `aux_3`, ordered CYP1A2, CYP2C9, CYP2D6, CYP3A4.

- The 6,146 challenge training compounds use existing recipient-excluded predictions:
  rows `outer_fold == fold` of the audited nested feature run
  `shared/20260912T052549384296Z-bac5fb3a`. Each is the mean of four pair models
  excluding that recipient's fold, including anchor recipients.
- The 750 test compounds use predictions from the single model refitted on all
  measured log2fc training data. This is not a fold-model average.
- `fold` and `role` are audit metadata, never additional model features.

`outputs/full.pt` stores the trained state, architecture and target scalers;
`membership.parquet` records every training compound. The full training history,
input snapshots, prior pair-model summaries and checkpoint provenance are retained.

The feature generator is `scripts/full_aux_log2fc.py`, launched with
`scripts/launch_full_aux_log2fc.py`. Independent verification is implemented in
`scripts/audit_full_aux_log2fc.py` and its launcher.

Independent audit `shared/20260912T074551656789Z-beecf03c` passed. All training
features are bitwise equal to their original recipient-excluded rows, training
membership equals all measured log2fc compound IDs, and reloading the saved model
reproduces all 750 test feature vectors with zero maximum absolute difference.
No test compound appears in training membership.
