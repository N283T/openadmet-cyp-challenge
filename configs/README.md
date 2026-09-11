# Experiment configurations

Place immutable run inputs under `regression/` or `classification/`; shared
feature-extraction configurations may live here. Use YAML or JSON, with explicit
seeds, model/hyperparameters, endpoint list, data snapshot IDs, feature version and
split ID. Pass the file with `cyp run --config FILE`; model code must load the saved
`CYP_CONFIG_PATH`. Do not infer changing defaults from a previous experiment.

Keep experiment-specific implementations under `src/cyp_challenge/` and small
entry scripts under `scripts/`. Add real configurations with the first baseline;
no placeholder model settings are selected by this infrastructure setup.

`regression/current-log2fc.json` is the explicit adopted auxiliary registry, not a
historical run configuration. `launch_parent_readout.py` reads it by default.
`regression/current-parent-readout.json` selects the corresponding five readout
runs used by `launch_tabicl_parent.py`. Pass `--revision` to allocate new immutable
experiment configurations. Historical comparisons can explicitly pass `--aux-log`
or `--readout-log`; they do not change the current registries.
