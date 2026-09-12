# Independent historical TabPFN inner reconstruction audit

Reconstruction run: `classification/20260912T132823500863Z-130ac6af`.
Independent audit: `classification/20260912T133039250120Z-68fe5a1c` — **passed**.

The audit reads frozen auxiliary targets, original UMAP assignments, ten auxiliary membership/metadata/history/checkpoint sets, all auxiliary prediction shards, frozen CheMeleon embeddings and Mordred tables, ten saved PCA/imputer objects and membership tables, and historical/new predictions as declared inputs. No additional model training is performed.

## Historical recipe agreement

Code review compares `classification_stack_inner.py` with the historical `descriptor_readout.py`, `stacking_readout.py`, `log2fc_pretrain.py` and `configs/screen/stack_mordred-fold0.json`.

The reconstructed inner learner preserves the historical feature/readout configuration: frozen2,048-dimensional CheMeleon embeddings, endpoint-training-only PCA256 without pre-PCA scaling, four nested out-of-sample log2fc predictions, and Mordred descriptors after training-only empty/constant removal and median imputation. The classifier remains TabPFN2.6 with eight estimators, seed42+outer fold, `fit_preprocessors` and the pinned classifier weights.

For outer fold `o`, inner fold `i=(o+1)%5`, and each remaining training fold `r`, that fold's log2fc features come from an encoder excluding `{o,i,r}`. Inner-validation log2fc predictions average the three corresponding encoders. Thus neither the assessed outer fold, inner fold, nor a readout-training row's own fold can contribute labels to its auxiliary feature model. Ten unique fold triples cover all thirty endpoint/outer/training-fold assignments.

The auxiliary fitting itself retains the historical seeded90/10 stopping split among eligible molecules and uses its best inner checkpoint. It does not add a full-eligible auxiliary refit that the historical pipeline did not use. Newer independently developed fine-tuning experiments use other refit policies; they are not silently substituted here.

## Artifact checks

| Check | Result |
|---|---|
| Ten auxiliary exclusions | Saved outer-holdout roles exactly equal the configured three-fold exclusions |
| Auxiliary train/stopping split | Exact reconstruction of seed42 permutation and90/10 partition |
| Four target means/stds | Recomputed from auxiliary training rows only; agreement within1e-7 |
| Auxiliary selected checkpoint | Epoch equals minimum recorded validation loss; checkpoint epoch/loss and SHA-256 match |
| Auxiliary predictions | Each source has6,896 unique molecule IDs and finite four-head predictions |
| Ten readout memberships | Exactly three training folds, designated inner fold, and excluded outer fold |
| PCA |256 components; fitted row count and mean exactly correspond to readout training data within1e-6 |
| Mordred filtering/imputation | Exact retained-column list and training medians; agreement within1e-12 |
| Inner predictions |5,081 canonical labels/IDs, correct endpoint/inner/outer fold keys, finite probabilities, recorded AUROC reproduced |
| Feature source assignments | All30 saved lineage entries match the required excluded triples |
| Historical predictions |5,081 outer OOF,7,500 fold-test and1,500 averaged-test probabilities exactly unchanged |

The completed inner outputs contain all1,497 D6 and3,584 3A4 observations, each used once as its designated inner validation observation. Their readout-training sample counts are859–963 D6 and2,111–2,193 3A4. Retained feature dimensions vary1,616–1,645 due to training-only descriptor filtering.

The emitted `inner-probabilities.parquet` is therefore suitable as the historical component's reconstructed inner source for the current fixed-pool ensemble comparison. Existing historical outer and test predictions remain separate, unchanged artifacts. This audit verifies their exact preservation and the inner reconstruction's source/mask/preprocessing provenance; it does not refit or replay TabPFN inference to manufacture another set of predictions.

The NumPy warnings about degrees of freedom arise when inspecting descriptors with no finite training values. Those columns are explicitly removed; the model input matrices are asserted finite by the training script. No audit assertion failed.

Files: `scripts/classification_stack_provenance_audit.py`, its launcher, `configs/classification-stack-inner/provenance-audit.json`, and the run's `outputs/audit.json`. Both audit scripts pass Ruff. No submission was made.
