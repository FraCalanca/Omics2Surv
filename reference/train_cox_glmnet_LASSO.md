# Train a Cox LASSO Model with Cross-Validation

Fits a Cox proportional hazards model using LASSO regularization (L1
penalty) with k-fold cross-validation and parallel processing support.

## Usage

``` r
train_cox_glmnet_LASSO(
  x_training,
  y_training,
  nfolds = 5,
  optimal_lambda = "minimum",
  type_measure = "C",
  penalty = NULL,
  system_type = "unix",
  nlambda = 200,
  seed = 123
)
```

## Arguments

- x_training:

  A matrix or data frame of predictors (e.g., omics data).

- y_training:

  A data frame containing survival outcomes, specifically
  `overall_survival` (time) and `status` (event).

- nfolds:

  Integer. Number of folds for cross-validation. Default is 5.

- optimal_lambda:

  Character. Choice of lambda for prediction: `"minimum"` (lambda.min)
  or `"1se"` (lambda.1se).

- type_measure:

  Character. Loss function for cross-validation. Default is `"C"`
  (C-index).

- penalty:

  Optional numeric vector for penalty factors (e.g., for Adaptive
  LASSO).

- system_type:

  Character. Operating system type, either `"unix"` or `"windows"`, to
  handle parallel backend registration.

- nlambda:

  Integer. The number of \\\lambda\\ values in the sequence. Default is
  200.

- seed:

  Integer. Seed for reproducibility of cross-validation folds.

## Value

A named list containing:

- `optimal_lambda`: The selection rule used.

- `cv_fit`: The fitted `cv.glmnet` object.

- `selected_lambda`: The numeric value of the chosen \\\lambda\\.

- `selected_features`: Character vector of predictors with non-zero
  coefficients.

- `pred`: Linear predictors for the training set.

- `cindex`: Training set Concordance Index.
