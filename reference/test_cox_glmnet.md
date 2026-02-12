# Evaluate Cox GLMNET Model on Test Data

Predicts risk scores and calculates the Concordance Index (C-index) for
a fitted Cox lasso model using independent test data.

## Usage

``` r
test_cox_glmnet(x_test, y_test, training_result)
```

## Arguments

- x_test:

  A matrix or data frame of predictor variables for the test set.

- y_test:

  A data frame or matrix containing survival information. Must include
  columns `overall_survival` (time) and `status` (event).

- training_result:

  A list containing the trained model objects:

  - `cv_fit`: A fitted `cv.glmnet` object.

  - `selected_lambda`: The penalty parameter (\\\lambda\\) to use for
    prediction.

## Value

A named list containing:

- `pred`: A matrix of predicted risk scores (linear predictors) for the
  test set.

- `cindex`: A numeric value representing the Harrell's Concordance
  Index.
