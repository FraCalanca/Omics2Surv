# Train Joint Multi-Omics Survival Models

Fits joint multi-omics survival models (`CoopCox`, `AFTCoop`, or
`blockForest`) on two or more omics layers from a training
`MultiAssayExperiment`. Survival outcome (`overall_survival` and
`status`) is extracted directly from `colData(mae_train)`. Results are
stored in the returned MAE.

## Usage

``` r
train_j_models(
  mae_train,
  omic_list,
  joint_method = c("CoopCox", "AFTCoop", "blockForest"),
  seed = 123,
  ...
)
```

## Arguments

- mae_train:

  A `MultiAssayExperiment` built with
  [`build_MAE`](https://fracalanca.github.io/Omics2Surv/reference/build_MAE.md).
  `colData` must contain `"overall_survival"` and `"status"`.

- omic_list:

  Character vector of assay names to use (e.g. `c("mRNA", "CNV")`).
  Order matters for `"AFTCoop"` (exactly 2 omics required).

- joint_method:

  Character string: `"CoopCox"` (2-3 omics), `"AFTCoop"` (exactly 2
  omics), or `"blockForest"` (any number).

- seed:

  Integer seed. Default is 123.

- ...:

  Additional arguments passed to the underlying train function:

  - `CoopCox`: `val_fraction`, `lambdaVector`, `folds`, `alpha`,
    `maxit`.

  - `AFTCoop`: `model`, `rho_values`, `ncore_max_rho`, `ncore_max_cv`.

  - `blockForest`: `threads`, `always.split.variables`,
    `n_top_features`.

## Value

`mae_train` with results stored at
`metadata(mae_train)$analysis$joint[[joint_method]]`, containing `train`
and `metrics` (resubstitution Harrell C-index).

## See also

[`predict_j_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_j_models.md),
[`train_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/train_cooplearning.md),
[`train_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/train_AFTCoop.md),
[`train_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/train_blockForest.md)
