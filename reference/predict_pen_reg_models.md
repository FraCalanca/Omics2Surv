# Predict from Trained Penalized Regression Models on a New MAE

Applies models trained by
[`train_pen_reg_models`](https://fracalanca.github.io/Omics2Surv/reference/train_pen_reg_models.md)
to a new `MultiAssayExperiment` (test set, validation cohort, or future
data). Survival outcome and clinical covariates are extracted directly
from `colData(mae_pred)` and `colData(mae_train)`.

## Usage

``` r
predict_pen_reg_models(
  mae_pred,
  mae_train,
  omic_list,
  integration = c("none", "early", "late"),
  model_method = c("cox_lasso", "cox_adaptive_lasso", "cox_elastic_net", "aft_pen"),
  clinical_vars = NULL
)
```

## Arguments

- mae_pred:

  A `MultiAssayExperiment` object containing the data to predict on,
  built with
  [`build_MAE`](https://fracalanca.github.io/Omics2Surv/reference/build_MAE.md)
  after aligning features via
  [`prepare_test_SE`](https://fracalanca.github.io/Omics2Surv/reference/prepare_test_SE.md).

- mae_train:

  A `MultiAssayExperiment` object containing fitted models in
  `metadata(mae_train)$analysis`, produced by
  [`train_pen_reg_models`](https://fracalanca.github.io/Omics2Surv/reference/train_pen_reg_models.md).

- omic_list:

  Character vector of assay names to predict on. Must match those used
  in
  [`train_pen_reg_models`](https://fracalanca.github.io/Omics2Surv/reference/train_pen_reg_models.md).

- integration:

  Character string specifying the integration strategy used during
  training: `"none"`, `"early"`, or `"late"`. Default is `"none"`.

- model_method:

  Character string specifying the model to retrieve: `"cox_lasso"`,
  `"cox_adaptive_lasso"`, `"cox_elastic_net"`, or `"aft_pen"`. Must
  match the `model_method` used in training. Default is `"cox_lasso"`.

- clinical_vars:

  Character vector of column names in `colData` to append as covariates.
  Must match those used during training. Default is `NULL`.

## Value

`mae_pred` with results stored in `metadata(mae_pred)$analysis`:

- `"none"` / `"early"`: per-assay results at
  `[[model_method]][[assay_name]]`, containing `test` (predictions) and
  `metrics` (Harrell + Uno C-index).

- `"late"`: per-omic results as above, plus ensemble results at
  `$late_integration[[model_method]][[tag]]` (where `tag` encodes the
  integrated omics), containing `test` and `metrics`. The `late_method`
  is retrieved automatically from `mae_train` metadata.

## See also

[`train_pen_reg_models`](https://fracalanca.github.io/Omics2Surv/reference/train_pen_reg_models.md),
[`prepare_test_SE`](https://fracalanca.github.io/Omics2Surv/reference/prepare_test_SE.md),
[`build_MAE`](https://fracalanca.github.io/Omics2Surv/reference/build_MAE.md)
