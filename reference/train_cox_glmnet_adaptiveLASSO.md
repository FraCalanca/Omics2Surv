# Train an Adaptive Cox LASSO Model

Fits a Cox proportional hazards model using Adaptive LASSO. This
two-step procedure uses Ridge regression coefficients to create penalty
factors, encouraging the model to select variables with stronger initial
signals.

## Usage

``` r
train_cox_glmnet_adaptiveLASSO(
  x_training,
  y_training,
  nfolds = 5,
  optimal_lambda = "minimum",
  type_measure = "C",
  unpenalized = NULL,
  system_type = "unix",
  nlambda = 150,
  seed = NULL
)
```

## Arguments

- x_training:

  A matrix or data frame of predictor variables.

- y_training:

  A data frame containing survival outcomes, specifically
  `overall_survival` (time) and `status` (event).

- nfolds:

  Integer. Number of folds for cross-validation in both steps.

- optimal_lambda:

  Character. Choice of lambda for prediction: `"minimum"` or `"1se"`.

- type_measure:

  Character. Loss function for cross-validation (e.g., `"C"`).

- unpenalized:

  Optional vector of indices or names for variables to be excluded from
  the penalty (penalty set to 0).

- system_type:

  Character. OS type: `"unix"` or `"windows"`.

- nlambda:

  Integer. Number of \\\lambda\\ values for the Ridge step sequence.

- seed:

  Optional integer for reproducibility.

## Value

A named list containing:

- `ridge_cv`: The preliminary Ridge cross-validation object.

- `cv_fit`: The final Adaptive LASSO cross-validation object.

- `selected_lambda`: The numeric \\\lambda\\ used for the final model.

- `selected_features`: Predictors with non-zero coefficients.

- `pred`: Linear predictors for the training set.

- `cindex`: Training set Concordance Index.

- `weights`: The penalty factors derived from the Ridge step.
