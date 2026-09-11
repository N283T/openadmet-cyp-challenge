# Matched-CYP Admetica fine-tuning

Fit run `20260911T151537020733Z-d635955c`; audit
`20260911T151844585530Z-6d02b666`. All 20 target/fold fine-tuned checkpoints,
last checkpoints, optimizer states, loss histories, memberships, adapted
embeddings and 20 downstream TabPFN states were saved.

## Protocol

Each CYP starts from its corresponding external Admetica 300-dimensional encoder.
A new 300->256->1 ReLU/dropout(.1) regression head predicts pIC50. Encoder and
head are updated together with AdamW, encoder LR 1e-5, head LR 1e-4, weight decay
1e-4, batch64, gradient norm clipping5. Loss is MSE of pIC50 standardized using
inner-training statistics only. No external auxiliary targets or other feature
columns are added; the pretrained encoder supplies the external information.

Current parent-retained outer folds are unchanged. Ten percent of the eligible
outer-training rows, randomly selected with seed42+fold, provide internal
validation; the remainder supplies gradient updates. Maximum60 epochs, patience10,
select minimum validation MSE and reload the best state. The direct prediction
is transformed back to pIC50 once. Some folds, including 2D6 and 3A4, select
near the60-epoch cap, so this initial schedule does not establish convergence.

A second readout fits TabPFN2.6/eight estimators on the adapted 300-dimensional
embeddings, using all eligible outer-training labels, including the inner
validation rows. This is allowed within the outer split; none of the outer
validation labels affects the encoder, checkpoint selection or TabPFN fitting.
The direct head and TabPFN have different effective gradient/readout sample sizes.

Audit checks exact molecule/endpoint/fold/label/bound agreement for 6,384 OOF
rows per variant. Metrics are saved and registered in PostgreSQL. No test
predictions or submission were produced.

## Results

| Model / fixed equal blend | All macro ST-RAE | Neighbor macro ST-RAE |
|---|---:|---:|
| Frozen Admetica emb + TabPFN | .743910 | .819000 |
| Fine-tuned direct head | .766816 | .891617 |
| Fine-tuned emb + TabPFN | .750941 | .848523 |
| FP+Mordred + pred_log2fc + frozen Admetica | **.661279** | .803584 |
| FP+Mordred + pred_log2fc + FT direct | .670587 | .830751 |
| FP+Mordred + pred_log2fc + FT emb | .661746 | .803868 |
| Above base pair + frozen Admetica + FT emb | .666842 | **.793020** |

Blends use equal component weights; no weight search. Replacing the frozen
component with the fine-tuned embedding component is approximately neutral in
the three-way ensemble. Keeping both trades worse all-OOF ST-RAE for better
neighbor ST-RAE.

| CYP | Frozen emb | Fine-tuned emb | Fine-tuned direct |
|---|---:|---:|---:|
| 1A2 | .800361 | .810056 | .830172 |
| 2C9 | .682876 | **.665736** | .684623 |
| 2D6 | .938661 | .963475 | .950349 |
| 3A4 | .553743 | .564498 | .602117 |

Fine-tuned embedding macro MAE .560771 versus frozen .557558; macro Spearman
.563083 versus .557726. Better rank correlation does not yield better ST-RAE
or MAE overall. The direct head has macro MAE .576601 and Spearman .544083.

Decision: keep the frozen Admetica component as the current baseline candidate.
Preserve the fine-tuned embedding model as an alternative with a neighbor-subset
ensemble benefit, but do not replace the current triple. This is one bounded
single-task pIC50 fine-tuning recipe, not a conclusion that all Admetica
fine-tuning or multitask/log2fc adaptation is ineffective. External pretraining
overlap remains unaudited, as in the frozen comparison.
