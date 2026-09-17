# Jeremy-medium source-family data collection

The public source families named for `chemprop_medium` are now represented
locally. This is a source acquisition milestone, not a claim that Jeremy's exact
47-head training matrix has been reconstructed. His public checkpoint and README
do not expose the original molecule manifest, task list, source revisions or
curation rules.

## New snapshots

### Full TDC ADME collection

Run **`20260916T081236163734Z-3ae85fcd`** snapshots all 27 datasets in the
current public `adme_dataset_names` registry at TDC revision
`c310c35f27e3f506411018ac43d97b8ba23ca652`. It also preserves that revision's
`metadata.py` and Harvard Dataverse V105 metadata for DOI
`10.7910/DVN/21LKWG`. The Dataverse metadata declares CC0-1.0.

All 29 HTTP requests succeeded on the first attempt. The 27 raw tabular files
contain **112,841 rows** in total. They cover absorption, distribution,
metabolism and excretion data, including the five Veith CYP inhibition panels,
three Carbon-Mangels substrate panels, clearance, half-life, PPB, solubility,
lipophilicity, permeability and related ADME measurements. No task selection,
thresholding, standardization, deduplication or challenge-overlap removal was
performed in this acquisition run.

Several TDC datasets overlap sources already acquired independently. In
particular, the AstraZeneca physicochemical/clearance data and Carbon-Mangels
substrate panels must not later be counted as independent copies merely because
they arrived through TDC.

### Polaris Novartis CYP3A4 surrogate panel

Run **`20260916T081236163734Z-aa990e05`** snapshots the public Polaris Novartis
CYP3A4 files plus the curation notebooks and README at `polaris-recipes`
revision `acab563e06f7d7cdc1aa9c927bfc9b6df7966219`. Polaris declares the dataset
CC-BY-4.0.

The raw surrogate training table has **16,373 rows**, the Zimmerlin test table
has **368 rows**, and Polaris's combined table has **16,740 rows**. The combined
data expose `log_kobs`, reversible-inhibition `pIC50_CYP3A4`, `pIC50_CYP2C9`
and `pIC50_CYP2D6` fields where available. The raw training parquet also
preserves 10 model columns per endpoint; those predictions must remain
distinguishable from measured labels during curation. The 16,740-row curated CSV
and all intermediate public parquets were retained instead of choosing one
representation prematurely.

## Source-family coverage after acquisition

| Public family named by Jeremy | Local status |
|---|---|
| ChEMBL | ChEMBL 36 mirror/extractions already available, including CYP and AZ ADME records |
| TDC | Full current public ADME registry now snapshotted (27 files) |
| Polaris | Primary-author Biogen six-endpoint panel already present; Novartis CYP3A4 panel now added |
| Novartis surrogate CYP | Polaris 16,740-row panel now present; the earlier EPA NVS five-CYP panel remains a distinct source |

Thus the practical data-availability blocker is removed at the **source-family**
level. Exact equality with Jeremy's corpus remains unknowable from public
artifacts. The next step is to standardize structures, exclude all challenge
connectivity blocks, identify duplicate publications across wrappers, and build
a compact endpoint matrix. The user explicitly deferred reproducing the unknown
47 task definitions, so no invented mapping is introduced here.
