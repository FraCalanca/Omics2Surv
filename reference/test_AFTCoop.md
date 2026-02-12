# Test AFT Cooperative Learning Model

Evaluates a trained AFT cooperative learning model on test data by
generating predictions and calculating concordance indices for each rho
value.

## Usage

``` r
test_AFTCoop(x_test1, x_test2, y_test, training_result)
```

## Arguments

- x_test1:

  A matrix or data frame of the first set of test predictors. Must have
  the same features (columns) as `x_training1` used in training.

- x_test2:

  A matrix or data frame of the second set of test predictors. Must have
  the same features (columns) as `x_training2` used in training.

- y_test:

  A data frame containing test survival outcomes with columns:

  - `overall_survival`: Numeric survival time

  - `status`: Numeric event indicator (1 = event, 0 = censored)

- training_result:

  A list object returned by
  [`train_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/train_AFTCoop.md)
  containing the trained model coefficients and parameters.

## Value

A list containing:

- pred:

  Named list of test set predictions for each rho value. Higher
  predicted values indicate higher risk (shorter predicted survival).

- cindex:

  Named list of concordance indices (C-index) on test data for each rho
  value. Values range from 0 to 1, where 0.5 indicates random
  predictions and 1.0 indicates perfect concordance.

## Details

This function applies the coefficients learned during training to new
test data. For each rho value in the trained model:

1.  Extracts the corresponding coefficient vector

2.  Generates risk predictions using `predict_aft_coop`

3.  Calculates the concordance index (C-index) comparing predictions to
    actual outcomes
