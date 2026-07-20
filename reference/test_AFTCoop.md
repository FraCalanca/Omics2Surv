# Test AFT Cooperative Learning Model

Applies a trained AFTCoop model to test data, generating predictions for
each rho value. The test matrices are centred using the training column
means (`mU`, `mZ`) stored by
[`train_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/train_AFTCoop.md),
so no test-set information leaks into the linear predictor.

## Usage

``` r
test_AFTCoop(x_test1, x_test2, training_result)
```

## Arguments

- x_test1:

  Numeric matrix (samples x features) for the first omics layer. Must
  match features used in training.

- x_test2:

  Numeric matrix (samples x features) for the second omics layer. Must
  match features used in training.

- training_result:

  List returned by
  [`train_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/train_AFTCoop.md).

## Value

A list with `pred`: named list of numeric prediction vectors, one per
rho value.

## See also

[`train_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/train_AFTCoop.md),
[`predict_j_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_j_models.md)
