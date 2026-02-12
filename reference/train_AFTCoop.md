# Train AFT Cooperative Learning Model

Trains an Accelerated Failure Time (AFT) model using cooperative
learning with two sets of predictors. The function fits a regularized
survival model across multiple correlation (rho) values and returns
coefficients, predictions, and performance metrics.

## Usage

``` r
train_AFTCoop(
  x_training1,
  x_training2,
  y_training,
  model = "weibull",
  rho_values = c(0.25, 0.5, 0.75),
  ncore_max_rho = 10,
  ncore_max_cv = 10,
  seed = 123
)
```

## Arguments

- x_training1:

  A matrix or data frame of the first set of training predictors (e.g.,
  first omic features).

- x_training2:

  A matrix or data frame of the second set of training predictors (e.g.,
  second omic features).

- y_training:

  A data frame containing survival outcomes with columns:

  - `overall_survival`: Numeric survival time

  - `status`: Numeric event indicator (1 = event, 0 = censored)

- model:

  Character string specifying the AFT distribution. Default is
  "weibull". Other options include "lognormal" and "loglogistic".

- rho_values:

  Numeric vector of correlation penalty parameters to evaluate. Default
  is `c(0.25, 0.5, 0.75)`. Controls the cooperation between the two
  predictor sets.

- ncore_max_rho:

  Integer specifying the maximum number of cores to use for parallel
  processing across rho values. Default is 10.

- ncore_max_cv:

  Integer specifying the maximum number of cores to use for
  cross-validation. Default is 10.

- seed:

  Integer seed for reproducibility. Default is 123.

## Value

A list containing:

- final_model_coefs:

  Data frame of coefficient estimates for each rho value. Rows
  correspond to features from both predictor sets, columns to rho
  values.

- sigma_est:

  Numeric estimated scale parameter (sigma) from the AFT model.

- rho_values:

  Numeric vector of rho values used in training (same as input).

- pred:

  Training set predictions for each rho value.

- cindex:

  Concordance indices (C-index) on training data for each rho value.

- selected_features:

  Named list of selected (non-zero) features for each rho value.

## Details

The function performs the following steps:

1.  Prepares survival data (time, event indicator, log-transformed time)

2.  Estimates the scale parameter (sigma) using
    [`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)

3.  Fits cooperative AFT models across specified rho values using
    `aft_coop`

4.  Generates predictions and calculates C-index for each rho value on
    training data

5.  Identifies selected features (non-zero coefficients) for each rho
    value
