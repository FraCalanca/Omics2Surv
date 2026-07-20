# Train a flexynesis model from a training MAE

MAE-level entry point: extracts the omics layers and survival outcome
from `mae_train`, trains flexynesis on the training cohort only (all
internal feature selection disabled), computes the resubstitution
Harrell C-index, and stores the fitted model, training predictions,
feature importance and metrics in
`metadata(mae_train)$analysis$flexynesis`.

## Usage

``` r
train_flexynesis_models(
  mae_train,
  omic_list,
  surv_time = "overall_survival",
  surv_event = "status",
  ...
)
```

## Arguments

- mae_train:

  A `MultiAssayExperiment`; `colData` must contain `"overall_survival"`
  and `"status"`.

- omic_list:

  Character vector of assay names to use.

- surv_time, surv_event:

  Survival column names. Defaults `"overall_survival"` / `"status"`.

- ...:

  Forwarded to
  [`train_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis.md)
  (e.g. `model_class`, `device`, `use_cv`, `threads`, `conda_env`,
  feature-selection overrides).

## Value

`mae_train` with results at `metadata(mae_train)$analysis$flexynesis`,
containing `train` (the
[`train_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis.md)
result), `metrics` (`harrell_train`) and `feature_importance`.

## See also

[`predict_flexynesis_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_flexynesis_models.md),
[`train_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis.md)
