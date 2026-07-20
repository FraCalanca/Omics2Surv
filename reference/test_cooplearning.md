# Test Cooperative Learning Model

Applies a trained CoopCox model to new data. Can be called directly or
via
[`predict_j_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_j_models.md).

## Usage

``` r
test_cooplearning(X_list_test, training_result)
```

## Arguments

- X_list_test:

  Named list of omics matrices (samples x features) for the test set.
  Must have the same omics layers and features, in the same order, as
  training.

- training_result:

  List returned by
  [`train_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/train_cooplearning.md).

## Value

A list with `pred`: named numeric vector of risk scores.

## Details

Omics layers and features are assumed to be in the same order as in
training.

## See also

[`train_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/train_cooplearning.md),
[`predict_j_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_j_models.md)
