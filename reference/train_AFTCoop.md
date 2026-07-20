# Train AFT Cooperative Learning Model

Trains an Accelerated Failure Time cooperative learning model on exactly
two omics layers. Called internally by
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md).

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

  Numeric matrix (samples x features) for the first omics layer.

- x_training2:

  Numeric matrix (samples x features) for the second omics layer.

- y_training:

  Data frame with `overall_survival` and `status`.

- model:

  Character string specifying the AFT distribution: `"weibull"`
  (default), `"lognormal"`, or `"loglogistic"`.

- rho_values:

  Numeric vector of cooperation penalty parameters to evaluate. Default
  is `c(0.25, 0.5, 0.75)`.

- ncore_max_rho:

  Integer maximum cores for parallelisation across rho values. Default
  is 10.

- ncore_max_cv:

  Integer maximum cores for cross-validation. Default is 10.

- seed:

  Integer seed. Default is 123.

## Value

A list with:

- `final_model_coefs`: Data frame of coefficients (features x rho
  values).

- `sigma_est`: Estimated AFT scale parameter.

- `rho_values`: Rho values used.

- `pred`: Named list of training predictions, one per rho.

- `selected_features`: Named list of selected features per rho.

- `mU`, `mZ`: Column means of the two training matrices, used to centre
  the test set consistently at prediction time.

## See also

[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md),
[`test_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/test_AFTCoop.md)
