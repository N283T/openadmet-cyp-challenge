# Jeremy's upstream medium encoder: frozen LGBM control

The actual public `chemprop_medium` encoder, used unchanged as a feature extractor, gives **0.651179 macro ST-RAE** with our fixed LightGBM readout. The neighboring-compound subset scores **0.739368**. Adding this model to the existing FP+Mordred / predicted-log2fc / Admetica trio with equal weights gives **0.642600** overall.

This is a diagnostic pretrained-model control. The author identifies the checkpoint as upstream public-data pretraining, and no evidence of challenge-label training was found. Its molecule manifest is unavailable, so the precise corpus cannot be independently audited; that absence is not evidence of leakage. It remains marked `formal_goal_eligible: false` pending the project's decision about eligibility of public pretrained controls.

## Which checkpoint was tested

Source repository revision: `d3a5549e0c98167435cbf8f5115c68a9a184ed7a`.

Checkpoint: `CYP_Challenge/checkpoints/chemprop_medium.pt`, SHA-256 `b8137a0d202669baa9747f21eac4ef78b000afe68a3f100d25c849c41ad8a3d5`.

The [source README](https://github.com/jeremycheminf/openadmet_scripts/blob/d3a5549e0c98167435cbf8f5115c68a9a184ed7a/CYP_Challenge/README.md) and `scripts/03_embed_frozen_encoders.py` identify this file as the sibling project's upstream public CYP/ADME pretrained checkpoint. The separate 17-head challenge fine-tuning procedure starts from this file; its adapted outputs were not used here. The README describes ChEMBL/TDC/Polaris and a public Novartis-surrogate CYP panel. This is broader pretraining than our scratch PubChem-only experiment.

Checkpoint metadata confirms width **600**, depth **4**, **NormAggregation(norm=100)**, and no batch normalization. Its original predictor has **47 heads**; task names are absent from its metadata and are not inferred here. The extraction preserves checkpoint graph transforms, which are identity transforms, and uses only message passing followed by the original aggregation.

## Recorded execution

| Stage | Run |
|---|---|
| Frozen extraction | `20260912T055819206227Z-c060d197` (shared) |
| Five-fold LightGBM readout | `20260912T055849013223Z-f16a5dfd` (regression) |
| Matched OOF / fixed-blend audit | `20260912T055910703420Z-770c2927` (regression) |

Loading used `weights_only=True` with the same explicitly allowlisted installed classes as the earlier metadata inspection; all encoder keys were loaded strictly. All 6,896 standardized challenge compounds produced finite 600-dimensional embeddings. The checkpoint, source README, extraction source, and compounds were declared snapshots. No challenge labels enter extraction and no encoder fine-tuning occurred.

Readout uses the current parent-retained five-fold CV and the established fixed LGBM recipe: L1 regression, 300 trees, 15 leaves, minimum child size 20, learning rate 0.05, L2 regularization 1. There is no additional descriptor, predicted-log2fc, or PCA input. The audit matched all **6,384** nonanchor compound/endpoint/fold keys and confidence intervals.

## Results

| Endpoint | ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| 1A2 | 0.728514 | 0.596155 | 0.533183 |
| 2C9 | 0.535216 | 0.409666 | 0.688735 |
| 2D6 | 0.858244 | 0.567901 | 0.444108 |
| 3A4 | 0.482742 | 0.508354 | 0.781702 |
| Macro | **0.651179** | **0.520519** | **0.611932** |

| Model / equal blend | All OOF | Neighbors |
|---|---:|---:|
| Actual medium + LGBM | 0.651179 | 0.739368 |
| Existing FP + log2fc + Admetica trio | 0.661279 | 0.803584 |
| Trio + actual medium | **0.642600** | 0.769872 |
| Trio + scratch PubChem600 | 0.658841 | 0.787963 |

The result supports testing a broader pretrained representation. It does not isolate whether its advantage comes from corpus breadth, particular source assays, normalization, depth, or optimization. This run contains no TabPFN comparison.

## Provenance timeline and code distinction

Recorded provenance check: shared run `20260912T060127278065Z-39e5a986`.

The complete local Git history shows one commit introducing the checkpoint: `cb09f313353e01e757ecf463e71f013c41c71d98`, authored and committed **2026-09-02 19:51:28 +01:00**. Its bytes at that first commit have the exact same SHA-256 as the checkpoint used here. Later README changes on September 2 and 4 did not replace the checkpoint. This establishes publication history, not the date of upstream training; publication was after the challenge began.

The public code makes the workflow distinction concrete:

- `scripts/03_embed_frozen_encoders.py` reads `checkpoints/chemprop_medium.pt` for frozen representations.
- `scripts/06_finetune_chemprop_multitask.py` uses that same file as `PRETRAINED`, constructs 17 challenge/external targets, and saves the fitted result separately to `results/chemprop_multitask_final_model.pt`.
- `src/cyp_submission/chemprop_transfer.py` loads the upstream model, retains its encoder and aggregation, and attaches a fresh target-specific head for the challenge adaptation.

Thus the tested 47-head checkpoint is the documented upstream input, not the public script's 17-head challenge-trained output. The author describes the upstream sources as public ChEMBL/TDC/Polaris and a Novartis-surrogate CYP panel. We did not find a separate explicit statement saying “challenge pIC50 excluded,” nor an upstream training timestamp or task/molecule manifest. No stronger assertion about individual upstream observations is made.
