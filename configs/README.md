# Experiment configurations

Place immutable run inputs under `regression/` or `classification/`; shared
feature-extraction configurations may live here. Use YAML or JSON, with explicit
seeds, model/hyperparameters, endpoint list, data snapshot IDs, feature version and
split ID. Pass the file with `cyp run --config FILE`; model code must load the saved
`CYP_CONFIG_PATH`. Do not infer changing defaults from a previous experiment.

Keep experiment-specific implementations under `src/cyp_challenge/` and small
entry scripts under `scripts/`. Add real configurations with the first baseline;
no placeholder model settings are selected by this infrastructure setup.
