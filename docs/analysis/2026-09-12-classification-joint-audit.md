# Independent audit of TDI joint fine-tuning

Completed model run `20260912T120135164747Z-1da37202` passed recorded audit
`20260912T120615145022Z-7e38be86`. The poor classification results do not come
from a detected held-out-label, target-scaling or export defect in this audit.

The audit reconstructs the two TDI heads and ten regression auxiliaries directly
from declared molecular measurements and auxiliary observations. For each of
the two recipes and five outer folds, it checks the membership masks over all
6,896 compounds in both phases. Inner training excludes the outer fold and
`(outer+1)%5`; outer refitting excludes only the outer fold. The same row mask
governs every auxiliary head. Test structures have no training membership.

Every saved auxiliary mean, population SD and observed training count exactly
matches a reconstruction using that phase's training rows. Classification heads
retain mean0/SD1. The audit replays the recorded inner-BCE improvement rule,
checks the selected epoch against `inner-best.pt`, and checks that `refit.pt`
and the actual refit history contain exactly that number of epochs. Source
inspection confirms each phase starts with the external-only encoder state and
a fresh head/optimizer; refit does not continue the inner checkpoint.

| Recipe | Selected epochs, outer folds0–4 |
|---|---|
| TDI only | 19,2,12,36,61 |
| TDI plus auxiliaries | 25,2,13,33,45 |

All20 endpoint/fold thresholds were independently recomputed from saved inner
probabilities using the source grid and tie policy. Outer IDs, folds, true
labels, applied thresholds and binary predictions agree. Each recipe has exactly
**5,081 canonical OOF predictions and1,500 mean test predictions**:1,497 D6 and
3,584 3A4 labels, with750 test compounds per endpoint. Five-fold test means
reconstruct with maximum absolute error zero. Exported test-label fields contain
no finite labels.

The audit checks saved output and training lineage, and reads checkpoint epoch
metadata; it does not rerun neural inference or independently replay optimization.
Its `fit-audit.csv` preserves per-fold epochs, thresholds and label counts, and
`summary.json` records the passed conditions. Scripts are
`scripts/classification_joint_audit.py` and
`scripts/classification_joint_audit_launch.py`; all inputs are snapshotted.
