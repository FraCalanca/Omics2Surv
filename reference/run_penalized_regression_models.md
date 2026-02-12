# Run Penalized Regression Models with Multi-Omics Integration

Applies penalized regression survival models (Lasso, adaptive Lasso,
elastic net, or AFT) using single omics, early integration (concatenated
features), or late integration (ensemble predictions).

## Usage

``` r
run_penalized_regression_models(
  omics_list,
  clinical_data,
  surv_list,
  integration = c("none", "early", "late"),
  model_method = c("cox_lasso", "cox_adaptive_lasso", "cox_elastic_net", "aft_pen"),
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
  omics layer), `"early"` (concatenate features before modeling), or
  `"late"` (model each layer separately, then ensemble).

- model_method:

  Character string specifying the penalized model: `"cox_lasso"`,
  `"cox_adaptive_lasso"`, `"cox_elastic_net"`, or `"aft_pen"`.

- train_fraction:

  Numeric proportion (0-1) of samples for training. Default is 0.7.

- seed:

  Integer seed for reproducibility. Default is 123.

- late_method:

  Character string for late integration: `"average"` or `"RRA"`. Only
  used when `integration = "late"`. Default is `"average"`.

- ...:

  Additional arguments passed to the penalized regression functions.

## Value

A MultiAssayExperiment object with results stored in
`metadata(mae)$analysis`:

- For `"none"` or `"early"`: Results at
  `metadata(mae)$analysis[[model_method]][[assay_name]]`

- For `"late"`: Results at
  `metadata(mae)$analysis$late_integration[[model_method]]`

Each contains `train` and `test` results with predictions and C-index.
