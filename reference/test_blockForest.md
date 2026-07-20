# Test Block Random Forest Model

Applies a trained block random forest model to new data. Can be called
directly or via
[`predict_j_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_j_models.md).

## Usage

``` r
test_blockForest(X_list, training_result, threads = 10)
```

## Arguments

- X_list:

  Named list of omics matrices (samples x features) for the test set.
  Names must match those used in training. Column names will be prefixed
  automatically as in training.

- training_result:

  List returned by
  [`train_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/train_blockForest.md).

- threads:

  Integer number of parallel threads. Default is 10.

## Value

A list with `pred`: numeric vector of risk scores.

## See also

[`train_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/train_blockForest.md),
[`predict_j_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_j_models.md)
