# Selected Features of a Penalised Model

Returns the features retained by a model that performs built-in feature
selection (e.g. the non-zero-coefficient predictors of `"cox_lasso"`,
`"cox_adaptive_lasso"`, `"cox_elastic_net"`), read from
`metadata(mae)$analysis[[model]][[modality]][[split]]$selected_features`.
The returned features can be passed to
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md).
For models *without* built-in selection (`"blockForest"`,
`"flexynesis"`), use
[`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)
instead.

## Usage

``` r
selected_features(mae, model, modality = NULL, split = "auto", rho = NULL)
```

## Arguments

- mae:

  A `MultiAssayExperiment` whose `metadata()$analysis` holds the fitted
  model.

- model:

  Character string naming the model (e.g. `"cox_lasso"`).

- modality:

  Character string naming the omics modality (e.g. an `"early_<omics>"`
  tag, `"mRNA"`, `"meth"`). Leave `NULL` for models stored without a
  modality level (joint methods), or to auto-resolve a model that has a
  single per-modality entry.

- split:

  Which split to read: `"auto"` (default; the single split present
  besides `"metrics"`, i.e. `"train"` in a training MAE) or an explicit
  `"train"` / `"test"`.

- rho:

  Numeric. For `model = "AFTCoop"` only: selects the `beta_rho_<rho>`
  feature set (e.g. `rho = 0.5`). If `NULL` (default) and only one rho
  was trained, it is resolved automatically.

## Value

A character vector of selected feature names.

## See also

[`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md),
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md)
