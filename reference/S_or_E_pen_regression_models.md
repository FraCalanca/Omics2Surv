# Single or Early Integration Penalized Regression Models

Trains and tests penalized regression survival models on a single omics
layer or early-integrated (concatenated) features, storing results in
the MultiAssayExperiment metadata.

## Usage

``` r
S_or_E_pen_regression_models(
  mae,
  assay_name,
  split,
  method = c("cox_lasso", "cox_adaptive_lasso", "cox_elastic_net", "aft_pen"),
  seed = 123,
  ...
)
```

## Arguments

- mae:

  A MultiAssayExperiment object.

- assay_name:

  Character string specifying the assay name to model (e.g., single
  omics layer or "early_integrated").

- split:

  A list containing:

  - `train_omic`: Training design matrix (samples × features)

  - `test_omic`: Test design matrix

  - `train_survival`: Training survival data frame with
    `overall_survival` and `status`

  - `test_survival`: Test survival data frame

- method:

  Character string specifying the penalized model: `"cox_lasso"`,
  `"cox_adaptive_lasso"`, `"cox_elastic_net"`, or `"aft_pen"`.

- seed:

  Integer seed for reproducibility. Default is 123.

- ...:

  Additional arguments passed to the training function (e.g., `alpha`,
  `nfolds`).

## Value

The input MultiAssayExperiment object with results stored at
`metadata(mae)$analysis[[method]][[assay_name]]` containing:

- `train`: Training results (model, predictions, C-index, selected
  features)

- `test`: Test results (predictions, C-index)
