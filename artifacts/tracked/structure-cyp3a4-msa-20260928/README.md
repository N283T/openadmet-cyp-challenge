# Full-length CYP3A4 alternative MSA

Recorded run: `shared/20260927T165502581117Z-7ca2b7c0`.
The A3M is in that run's `outputs/P08684.a3m` (8,300,528 bytes), with
SHA-256 `809d4744a32a1291b14e2940d88662e7e68389a6872353a9fb707737ca219071`.
The external A3M payload remains local; this directory tracks its compact
provenance and the source manifest row.

Source: [LiteFold/human-proteome-wide-msa](https://huggingface.co/datasets/LiteFold/human-proteome-wide-msa),
revision `632cfc87f567e429570897c3addd80665a697f3e`, file
[`a3m/shard-002/P08684.a3m`](https://huggingface.co/datasets/LiteFold/human-proteome-wide-msa/blob/632cfc87f567e429570897c3addd80665a697f3e/a3m/shard-002/P08684.a3m).
The dataset card states that it was generated with ColabFold/MMseqs2 GPU using
UniRef30 2023-02 and the ColabFold environmental database 2021-08, under
CC BY 4.0.

The first A3M sequence exactly matches the official P08684 503-residue
sequence. The file contains 14,300 aligned sequences, and every A3M row spans
503 query columns after removing lowercase insertions. This is a separately
generated MSA, **not** the AFDB v6 MSA used for the reference prediction.
The AFDB API URL still returned HTTP 403 when checked on 2026-09-28.
