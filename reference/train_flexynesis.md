# Train a flexynesis model on the training cohort only

Exports the training cohort, runs flexynesis in training mode with all
feature selection disabled, and persists the model and artifacts for
later inference. The test cohort is never seen here. flexynesis still
requires a `test/` folder to start, so a copy of the training cohort is
written as a harmless placeholder (identical features =\> the harmonize
step is a no-op and nothing from a real test set can leak into
training).

## Usage

``` r
train_flexynesis(
  X_train_list,
  y_train,
  omic_list,
  surv_time = "overall_survival",
  surv_event = "status",
  tmpdir = tempfile("flexynesis_"),
  conda_env = "flexynesisenv",
  conda_path = NULL,
  check_install = TRUE,
  safetensors = FALSE,
  ...
)
```

## Arguments

- X_train_list:

  Named list of training omics matrices (samples x features).

- y_train:

  Data frame of training survival data; row names are sample IDs.

- omic_list:

  Character vector of omics layer names to use.

- surv_time, surv_event:

  Survival column names. Defaults `"overall_survival"` / `"status"`.

- tmpdir:

  Working directory; a fresh per-call subdirectory by default.

- conda_env, conda_path:

  Conda environment / binary.

- check_install:

  Logical; verify flexynesis importable first (default TRUE).

- safetensors:

  Logical; passed through to
  [`run_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/run_flexynesis.md).

- ...:

  Additional arguments forwarded to
  [`run_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/run_flexynesis.md)
  (e.g. `model_class`, `device`, `use_cv`, `threads`, or the
  feature-selection overrides).

## Value

A list with `model` (paths to the saved model/artifacts and the working
dirs), `pred` (training-split predictions), `feature_importance`
(stacked importance table, or `NULL`), `selected_features` (alias of
`feature_importance`), and `params`.

## See also

[`test_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/test_flexynesis.md),
[`run_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/run_flexynesis.md)
