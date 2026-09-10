#!/usr/bin/env bash
set -euo pipefail
for fold in 0 1 2 3 4; do
  pixi run cyp run --track shared \
    --purpose "CheMeleon-initialized log2fc pretraining; outer UMAP fold ${fold} fully excluded; inner-only scaling and best-checkpoint selection" \
    --config "configs/pretraining/chemeleon-log2fc-fold${fold}.json" \
    --input artifacts/runs/shared/20260910T100337258849Z-927eac1c/outputs/targets.parquet \
    --input artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet \
    --input artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet \
    --input artifacts/runs/shared/20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt \
    --input artifacts/runs/shared/20260910T012612848389Z-bfeed434/outputs/provenance.json \
    -- python -m cyp_challenge.log2fc_pretrain
done
