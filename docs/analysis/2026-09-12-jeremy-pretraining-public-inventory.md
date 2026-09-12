# Jeremy's upstream 47-task pretraining: bounded public-source inventory

The **47-task definition, training compound manifest and dataset assembly code
were not recovered** from the checked public sources. This is a bounded negative
finding, not evidence that the project never existed or that its data were private.
The two usable upstream checkpoint files remain publicly available.

## Scope and recorded evidence

- Initial snapshot run: **`20260912T072002498987Z-82e593f0`**. The unauthenticated
  [author public repository listing](https://api.github.com/users/jeremycheminf/repos?per_page=100&page=1)
  returned 27 repositories in one complete page. No sibling pretraining repository
  appears by name. Candidate public repository endpoints `adme_pretrain`,
  `chemprop_pretraining`, `chemprop_pretrain`, and `adme-pretrain` each returned
  HTTP 404. Public 404 responses cannot distinguish nonexistent from private repos.
- Detail snapshot run: **`20260912T072115632580Z-51dac613`**. Inspected complete
  recursive source trees of `openadmet_scripts`, the author's blog,
  `colab_notebooks`, `streamlit_apps`, and the author's `chemprop` fork. No tree
  response was truncated. The blog contains a welcome post; the notebook repository
  contains SureChEMBL similarity and ChEMBL feature-generation notebooks; the app
  repository contains correlation and SMIRKS tools. None exposes a named upstream
  ADME pretraining project in its source tree.
- [CYP directory history](https://api.github.com/repos/jeremycheminf/openadmet_scripts/commits?path=CYP_Challenge&per_page=100)
  contains five commits. All five historical README files and complete repository
  trees were fetched. The first CYP commit is
  `cb09f313353e01e757ecf463e71f013c41c71d98`; the current revision is
  `d3a5549e0c98167435cbf8f5115c68a9a184ed7a`. The public branch listing contains
  only `main`. The initial kit already describes the sibling pretraining project
  without providing its task list or source directory. Subsequent commits change
  attribution, submission handling and calibration, not upstream data assembly.
- Every fetched HTTP body, response status, URL and SHA256 is preserved in the
  corresponding run. Scripts are `jeremy_pretrain_inventory_fetch.py` and
  `jeremy_pretrain_inventory_detail.py`; no training or large data download occurred.

## What is actually available

The [CYP README](https://github.com/jeremycheminf/openadmet_scripts/blob/d3a5549e0c98167435cbf8f5115c68a9a184ed7a/CYP_Challenge/README.md)
attributes pretraining data broadly to public ChEMBL, TDC, Polaris and a public
Novartis-surrogate CYP panel. These are source families, not a versioned dataset
manifest. Its references provide ChEMBL/Capricho curation background and an
OpenADMET public-data talk, without mapping 47 output columns to datasets.

[Script03](https://github.com/jeremycheminf/openadmet_scripts/blob/d3a5549e0c98167435cbf8f5115c68a9a184ed7a/CYP_Challenge/scripts/03_embed_frozen_encoders.py)
identifies two frozen upstream encoders from the same sibling pretraining project:
medium D-MPNN and a CheMeleon-initialized model. Their public checkpoint tensors
verify 47 heads, but earlier recorded checkpoint inspections found no named output
columns or training molecule manifest. The included `data/external/cyp_chembl.csv`
is challenge auxiliary CYP data; it must not be mistaken for the full 47-task
pretraining table. The challenge's 17-head fine-tuning script is similarly a
separate downstream training stage.

## Consequence for our comparison

We can compare the actual frozen upstream representations under our fixed CV,
which has been done. We **cannot attribute medium's advantage over our 43,496-
compound corpus to a larger upstream sample count**, because that count and the
composition are unknown. Endpoint breadth, assay pooling, sampling, training
schedule and architectural differences are all unresolved explanations.

The concrete reproducible assets are the two upstream checkpoints, their frozen
extraction code, and the downstream CYP auxiliary data/code. Reconstructing the
original 47-task corpus requires an additional public release of the task schema
and assembly scripts or molecule/source manifest. No such release was found in
this bounded author-repository and history inspection; no message was sent to the
author.
