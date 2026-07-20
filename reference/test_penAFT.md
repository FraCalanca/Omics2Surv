# Evaluate a Penalized AFT Model on Test Data

Computes linear predictors for a fitted penalized AFT model on
independent test data, using the cross-validated fit produced by
[`train_penAFT`](https://fracalanca.github.io/Omics2Surv/reference/train_penAFT.md).

## Usage

``` r
test_penAFT(x_test, training_result)
```

## Arguments

- x_test:

  A matrix or data frame of predictors for the test set (samples in
  rows, same features as the training set).

- training_result:

  A list returned by
  [`train_penAFT`](https://fracalanca.github.io/Omics2Surv/reference/train_penAFT.md).
  Must contain `best_cv_fit`, the `penAFT.cv` object used for
  prediction.

## Value

A named list with:

- `pred`: Linear predictors for the test set, named by the row names of
  `x_test`.

## See also

[`train_penAFT`](https://fracalanca.github.io/Omics2Surv/reference/train_penAFT.md)
