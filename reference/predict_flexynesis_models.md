# Predict from a trained flexynesis model on a new MAE

Retrieves the model fitted by
[`train_flexynesis_models`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis_models.md)
from `mae_train`, runs flexynesis inference on `mae_pred`, and computes
the Harrell and Uno C-indices (the latter using the training cohort for
IPCW, exactly as the joint methods do). Results are stored in
`metadata(mae_pred)$analysis$flexynesis`.

## Usage

``` r
predict_flexynesis_models(mae_pred, mae_train, omic_list, ...)
```

## Arguments

- mae_pred:

  A `MultiAssayExperiment` to predict on.

- mae_train:

  A `MultiAssayExperiment` carrying the fitted model in
  `metadata(mae_train)$analysis$flexynesis`.

- omic_list:

  Character vector of assay names; must match training.

- ...:

  Forwarded to
  [`test_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/test_flexynesis.md).

## Value

`mae_pred` with results at `metadata(mae_pred)$analysis$flexynesis`,
containing `test` and `metrics` (`harrell_test` + `uno_test`).

## See also

[`train_flexynesis_models`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis_models.md),
[`test_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/test_flexynesis.md)
