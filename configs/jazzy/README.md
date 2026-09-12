# Jazzy experiments

Install the isolated, locked extractor environment with `pixi install -e jazzy`.
Jazzy 0.1.4 pins RDKit 2024.3.1; the default model environment stays unchanged.
Run extraction through the default experiment recorder, using the absolute
interpreter path because recorded commands execute in their source snapshot:

```sh
pixi run cyp run --track shared \
  --purpose 'Generate standard Jazzy ETDG 3-seed descriptors' \
  --config configs/jazzy/features_etdg_3seeds.json \
  --input artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet \
  --input pixi.lock \
  -- "$PWD/.pixi/envs/jazzy/bin/python" scripts/jazzy_features.py
```

Extraction uses standardized SMILES and the package's standard ETDG implementation
with seeds 11, 42, and 89, without force-field minimization or 2D fallback.
Descriptors are `sdc`, `sdx`, `sa`, `dga`, `dgp`, and `dgtot`. All seed results,
failures, molecule IDs, extraction timing, and package versions are retained.
Mean and population standard deviation use successful seeds; zero successes
produce missing values. These three seeds measure local seed sensitivity rather
than a thermodynamic conformer ensemble.

The public API fixes Kallisto total charge to zero. Net-charged molecules and
disconnected structures therefore receive missing features, without changing
their source SMILES. The downstream comparison retains those molecules and fits
imputation only on outer-training rows.

The launcher `scripts/launch_jazzy_readout.py` records paired readouts and an audit;
each invocation requires a unique `--tag` to preserve completed configurations.
Its feature argument must identify a completed extraction artifact, with six
`--mean-columns` and six `--std-columns` as declared in the extraction output.
The default backend is LightGBM; `--backend tabpfn` runs the matched TabPFN 2.6
comparison and declares its checkpoint as an input. Saved configurations contain
the exact column order, baseline reference, and permitted reproduction tolerance.
