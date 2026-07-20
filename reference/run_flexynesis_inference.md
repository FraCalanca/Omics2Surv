# Run the flexynesis CLI (inference mode)

Predicts on a test-only cohort using a model + artifacts produced by
[`run_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/run_flexynesis.md).
Inference forces CPU inside flexynesis and reuses the exact features and
scalers stored in the artifacts, so the feature selection flags do not
apply here.

## Usage

``` r
run_flexynesis_inference(
  pretrained_model,
  artifacts,
  data_path_test,
  outdir,
  surv_time = "overall_survival",
  surv_event = "status",
  prefix = "job",
  safetensors = FALSE,
  extra_args = character(0),
  conda_env = "flexynesisenv",
  conda_path = NULL
)
```

## Arguments

- pretrained_model:

  Path to `.final_model.pth` / `.safetensors`.

- artifacts:

  Path to `.artifacts.joblib` / `.artifacts.json`.

- data_path_test:

  FLAT folder with `clin.csv` + omics CSVs (see
  [`export_flexynesis_data`](https://fracalanca.github.io/Omics2Surv/reference/export_flexynesis_data.md)
  with `layout = "flat"`).

- outdir:

  Output folder for predictions.

- surv_time, surv_event:

  Survival column names. These MUST be passed: in flexynesis the
  predicted-labels block is gated by
  `if any([target_variables, surv_event_var])`, so without
  `--surv_event_var` inference stops after writing embeddings and never
  produces `predicted_labels.csv`.

- prefix:

  Job prefix. Default `"job"`.

- safetensors:

  Logical; set if the model was saved as SafeTensors.

- extra_args:

  Optional additional raw CLI tokens.

- conda_env, conda_path:

  Conda environment / binary.

## Value

Invisibly the path to `<prefix>.predicted_labels.csv`.
