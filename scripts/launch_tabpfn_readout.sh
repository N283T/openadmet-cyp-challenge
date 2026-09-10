#!/usr/bin/env bash
set -euo pipefail
parents=(20260910T100719601351Z-d8f745ca 20260910T100846602294Z-f7a00ae9 20260910T101011080180Z-d3607e4f 20260910T101136302254Z-e2a36516 20260910T101301049877Z-5d886b77)
for fold in 0 1 2 3 4; do
  parent="artifacts/runs/shared/${parents[$fold]}/outputs"
  args=()
  for file in "$parent"/embeddings-*.parquet; do args+=(--input "$file"); done
  pixi run cyp run --track shared --purpose "TabPFN-2.6 on fold-matched log2fc-pretrained embeddings only; UMAP fold ${fold}" \
    --config "configs/readout/tabpfn-embedding-fold${fold}.json" \
    --input artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet \
    --input artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet \
    --input artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet \
    --input "$parent/summary.json" --input "$parent/membership.parquet" \
    --input /home/nagaet/.cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt \
    --input /home/nagaet/.cache/tabpfn/tabpfn-v2.6-classifier-v2.6_default.ckpt \
    "${args[@]}" -- python scripts/tabpfn_readout.py
done
