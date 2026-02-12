# Run Model Comparison Pipeline for Multi-Omics Survival Analysis

Executes a comprehensive comparison of multiple survival models across
different integration strategies (none, early, late, or joint),
returning performance metrics for each model.

## Usage

``` r
run_model_comparison_pipeline(
  omics_list,
  clinical_data,
  surv_list,
  integration = c("none", "early", "late", "joint"),
  model_methods,
  train_fraction = 0.7,
  seed = 123,
  late_method = c("average", "RRA"),
  ...
)
```

## Arguments

- omics_list:

  Named list of SummarizedExperiment objects, each representing an omics
  layer.

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns).

- surv_list:

  Character vector of survival variable names (e.g.,
  `c("overall_survival", "status")`).

- integration:

  Character string specifying integration strategy: `"none"` (single
  omics), `"early"` (concatenated features), `"late"` (ensemble
  predictions), or `"joint"` (cooperative modeling).

- model_methods:

  Character vector of model names to compare. Options depend on
  `integration`: for "none"/"early", use penalized regression models;
  for "late", use same; for "joint", use
  `c("CoopCox", "AFTCoop", "blockForest", "flexynesis")`.

- train_fraction:

  Numeric proportion (0-1) of samples for training. Default is 0.7.

- seed:

  Integer seed for reproducibility. Default is 123.

- late_method:

  Character string for late integration: `"average"` or `"RRA"`. Only
  used when `integration = "late"`. Default is `"average"`.

- ...:

  Additional arguments passed to model-specific functions.

## Value

A list containing:

- `comparison`: Data frame with columns `model`, `train_cindex`, and
  `test_cindex` for each model tested

- `per_model_results`: Named list of detailed results for each model,
  including reports from `metadata(mae)$analysis`
