# Run the flexynesis CLI (training mode)

Builds and executes a flexynesis training command via `conda run`, with
all internal feature selection disabled by default so pre-filtered data
passes through untouched. Output is streamed live; on completion the
function checks the log for dropped samples and verifies the expected
output files exist.

## Usage

``` r
run_flexynesis(
  data_path,
  omic_list,
  surv_time,
  surv_event,
  outdir,
  model_class = "DirectPred",
  prefix = "job",
  threads = 10,
  device = "cpu",
  use_cv = TRUE,
  variance_threshold = 0,
  correlation_threshold = 1,
  features_top_percentile = 100,
  features_min = 0,
  subsample = 0,
  finetuning_samples = 0,
  feature_importance_method = c("IntegratedGradients", "GradientShap", "Both"),
  safetensors = FALSE,
  extra_args = character(0),
  conda_env = "flexynesisenv",
  conda_path = NULL
)
```

## Arguments

- data_path:

  Folder holding `train/` and `test/` (see
  [`export_flexynesis_data`](https://fracalanca.github.io/Omics2Surv/reference/export_flexynesis_data.md)
  with `layout = "train_test"`).

- omic_list:

  Character vector of omics layer names; becomes `--data_types` and must
  match the CSV file names.

- surv_time, surv_event:

  Column names in `clin.csv` for survival time and event indicator.

- outdir:

  Output folder for flexynesis results.

- model_class:

  flexynesis model class. Default `"DirectPred"`. With
  `surv_time`/`surv_event` set, the DL models (DirectPred,
  supervised_vae, MultiTripletNetwork, CrossModalPred, GNN) attach an
  MLP head trained with Cox-PH loss, so they predict survival risk
  scores directly.

- prefix:

  Job prefix for output files. Default `"job"`.

- threads:

  Number of CPU threads. Default 10.

- device:

  `"cpu"`, `"mps"`, `"cuda"` or `"auto"`.

- use_cv:

  Logical; request 5-fold cross-validation.

- variance_threshold, correlation_threshold, features_top_percentile,
  features_min:

  Feature-selection controls. Defaults (0, 1, 100, 0) keep every
  feature, i.e. no selection. `features_top_percentile = 100` (not 0) is
  used on purpose: 0 makes older flexynesis skip the selection log and
  crash with `KeyError: 'select_features'`, whereas 100 keeps all
  features via the Laplacian early-return. Override only if you want
  flexynesis to filter.

- subsample:

  Training downsample size; 0 (default) keeps every patient.

- finetuning_samples:

  Test samples used for fine-tuning; 0 (default) off.

- feature_importance_method:

  Captum explainer for marker discovery: `"IntegratedGradients"`
  (default), `"GradientShap"` or `"Both"`. Marker finding is left
  enabled so the importance CSV(s) are written.

- safetensors:

  Logical; save as SafeTensors + JSON instead of .pth + .joblib.

- extra_args:

  Optional character vector of additional raw CLI tokens.

- conda_env, conda_path:

  Conda environment / binary.

## Value

Invisibly a list of output paths: `model`, `artifacts`,
`predicted_labels`, `feature_importance` (character vector of the
importance CSV(s) actually written).
