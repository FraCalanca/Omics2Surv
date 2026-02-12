# Train a Penalized Accelerated Failure Time (AFT) Model

Fits a penalized AFT model using the `penAFT` package. The function
optimizes the `alpha` parameter through cross-validation to balance
regularization types, supporting both standard Elastic Net and Sparse
Group LASSO.

## Usage

``` r
train_penAFT(
  x_training,
  y_training,
  nfolds = 5,
  penalty = "EN",
  alpha_range = seq(0.1, 0.9, by = 0.1),
  groups = NULL,
  nlambda = 100,
  seed = 123
)
```

## Arguments

- x_training:

  A matrix or data frame of predictors (samples in rows).

- y_training:

  A data frame containing survival data with `overall_survival` and
  `status` columns.

- nfolds:

  Integer. Number of folds for cross-validation. Default is 5.

- penalty:

  Character. The penalty type to apply: `"EN"` (Elastic Net) or `"SG"`
  (Sparse Group LASSO).

- alpha_range:

  Numeric vector. The range of `alpha` values to test (mixture between
  L1 and L2, or between Group and individual LASSO).

- groups:

  A vector of integers indicating group membership for predictors.
  Required if `penalty = "SG"`.

- nlambda:

  Integer. Number of \\\lambda\\ values in the tuning grid.

- seed:

  Integer or NULL. Seed for reproducibility.

## Value

A named list containing:

- `best_cv_fit`: The `penAFT.cv` object for the optimal `alpha`.

- `best_alpha`: The selected `alpha` value.

- `selected_lambda`: The value of \\\lambda\\ that minimized CV error.

- `selected_features`: Names of features with non-zero coefficients.

- `coefficients`: Vector of estimated coefficients.

- `pred`: Linear predictors for the training set.

- `cindex`: Training set Concordance Index.

- `cv_results`: A list of results for all tested `alpha` values.
