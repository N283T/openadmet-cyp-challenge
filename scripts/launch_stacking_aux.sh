#!/usr/bin/env bash
set -euo pipefail
for config in configs/stacking/aux-*.json; do
  pixi run cyp run --track shared --purpose "Nested log2fc cross-fitting: exclude UMAP pair ${config}" \
    --config "$config" \
    --input artifacts/runs/shared/20260910T100337258849Z-927eac1c/outputs/targets.parquet \
    --input artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet \
    --input artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet \
    --input artifacts/runs/shared/20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt \
    --input artifacts/runs/shared/20260910T012612848389Z-bfeed434/outputs/provenance.json \
    -- python -m cyp_challenge.log2fc_pretrain
done
