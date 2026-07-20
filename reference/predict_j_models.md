# Predict from Trained Joint Multi-Omics Models

Applies joint models trained by
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md)
to a new MAE. The underlying test functions
([`test_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/test_cooplearning.md),
[`test_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/test_AFTCoop.md),
[`test_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/test_blockForest.md))
can always be called directly.

## Usage

``` r
predict_j_models(
  mae_pred,
  mae_train,
  omic_list,
  joint_method = c("CoopCox", "AFTCoop", "blockForest"),
  ...
)
```

## Arguments

- mae_pred:

  A `MultiAssayExperiment` to predict on.

- mae_train:

  A `MultiAssayExperiment` with fitted models in
  `metadata(mae_train)$analysis$joint`.

- omic_list:

  Character vector of assay names. Must match those used in
  [`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md).

- joint_method:

  Character string: `"CoopCox"`, `"AFTCoop"`, or `"blockForest"`.

- ...:

  Additional arguments passed to the test function (e.g. `threads` for
  blockForest).

## Value

`mae_pred` with results stored at
`metadata(mae_pred)$analysis$joint[[joint_method]]`, containing `test`
and `metrics` (Harrell + Uno C-index).

## See also

[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md),
[`test_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/test_cooplearning.md),
[`test_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/test_AFTCoop.md),
[`test_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/test_blockForest.md)
