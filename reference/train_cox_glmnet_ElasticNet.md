# Train a Cox Elastic Net Model with Alpha Optimization

Fits a Cox proportional hazards model using Elastic Net regularization.
The function iterates through a range of `alpha` values to identify the
optimal mixture of L1 and L2 penalties based on cross-validated
performance.

## Usage

``` r
train_cox_glmnet_ElasticNet(
  x_training,
  y_training,
  nfolds = 5,
  optimal_lambda = "minimum",
  type_measure = "C",
  alpha_range = seq(0.1, 0.9, by = 0.1),
  penalty = NULL,
  system_type = "unix",
  nlambda = 200,
  seed = 123
)
```

## Arguments

- x_training:

  A matrix or data frame of predictor variables.

- y_training:

  A data frame containing survival outcomes with columns
  `overall_survival` and `status`.

- nfolds:

  Integer. Number of folds for cross-validation. Default is 5.

- optimal_lambda:

  Character. Rule for lambda selection: `"minimum"` (lambda.min) or
  `"1se"` (lambda.1se).

- type_measure:

  Character. Performance metric for CV (e.g., `"C"` for C-index or
  `"deviance"`).

- alpha_range:

  Numeric vector. The sequence of `alpha` values (0 to 1) to test.
  Default is 0.1 to 0.9.

- penalty:

  Optional numeric vector of penalty factors for individual features.

- system_type:

  Character. Operating system: `"unix"` or `"windows"`.

- nlambda:

  Integer. Number of \\\lambda\\ values in each `glmnet` path. Default
  is 200.

- seed:

  Integer. Seed for reproducibility.

## Value

A named list containing:

- `cv_fit`: The `cv.glmnet` object for the best performing `alpha`.

- `best_alpha`: The alpha value that yielded the best CV result.

- `selected_lambda`: The numeric \\\lambda\\ used for the final model.

- `selected_features`: Predictors with non-zero coefficients at the
  chosen alpha/lambda.

- `pred`: Linear predictors for the training set.

- `cindex`: Training set Concordance Index.
