# Stratify Patients into Risk Groups from Training and Test MAEs

Single entry point that (1) extracts a model's per-patient risk scores
from the training `MultiAssayExperiment`, (2) derives a cut-point on
those training scores, and (3) extracts and stratifies the test
predictions using that cut-point. Survival outcome (`overall_survival`
and `status`) is read directly from `colData(mae)`; extraction and
coercion are handled internally. Predictions are read from
`metadata(mae)$analysis[[model]][[modality]][[split]]$pred`.

## Usage

``` r
patient_stratification(
  mae_train,
  mae_test = NULL,
  model,
  modality = NULL,
  container = NULL,
  rho = NULL,
  method = c("median", "mediana", "user", "survcutpoint"),
  threshold = NULL,
  plot_scores = FALSE,
  minprop = 0.1,
  inclusive = FALSE,
  colors = .risk_colors(),
  time_col = "overall_survival",
  status_col = "status"
)
```

## Arguments

- mae_train:

  A `MultiAssayExperiment` holding the training split; the threshold and
  the optional score histogram are derived from it. `colData` must
  contain `"overall_survival"` and `"status"`.

- mae_test:

  A `MultiAssayExperiment` holding the test split, to be stratified with
  the training threshold. If `NULL`, the threshold is applied back to
  the training set.

- model:

  Character string naming the model in `metadata()$analysis` (e.g.
  `"cox_elastic_net"`).

- modality:

  Character string naming the omics modality (e.g. an `"early_<omics>"`
  tag, `"mRNA"`, `"meth"`). Leave `NULL` for models stored without a
  modality level (joint methods such as `"blockForest"`, and
  `"flexynesis"`), or to auto-resolve a model that has a single
  per-modality entry (e.g. one early-integration tag, or the single
  ensemble under `container = "late_integration"`).

- container:

  Optional container to disambiguate models whose names are reused
  across integration strategies. Use `"late_integration"` to read the
  consensus prediction of a late-integration ensemble (e.g.
  `model = "cox_lasso"`, `container = "late_integration"`); the sample
  IDs are taken from its `merged_risks` table.

- rho:

  Numeric. For `model = "AFTCoop"` only: selects the `beta_rho_<rho>`
  variant (e.g. `rho = 0.5`). If `NULL` (default) and only one rho was
  trained, it is resolved automatically; if multiple rho values are
  present, an error lists the available keys.

- method:

  Threshold method:

  - `"median"` (alias `"mediana"`): median of the training scores.

  - `"user"`: user-supplied value, passed via `threshold`.

  - `"survcutpoint"`: optimal cut-point via
    [`surv_cutpoint`](https://rdrr.io/pkg/survminer/man/surv_cutpoint.html).

- threshold:

  Numeric cut-point. Required when `method = "user"`.

- plot_scores:

  Logical. If `TRUE`, draws and returns a histogram of the training
  scores coloured by risk group, with the threshold marked. Default is
  `FALSE`.

- minprop:

  Minimum group proportion for
  [`surv_cutpoint`](https://rdrr.io/pkg/survminer/man/surv_cutpoint.html).
  Default is 0.1.

- inclusive:

  Logical. If `TRUE`, scores equal to the threshold are assigned to
  `"High risk"`. Default is `FALSE`.

- colors:

  Named character vector mapping `"High risk"` and `"Low risk"` to
  colours. Default is
  `c("High risk" = "red", "Low risk" = "deepskyblue")`.

- time_col, status_col:

  Names of the survival time and event columns in `colData`. Default to
  `"overall_survival"` and `"status"`.

## Value

A `list` with:

- `strata`: the stratified test `data.frame` (pass to
  [`plot_km`](https://fracalanca.github.io/Omics2Surv/reference/plot_km.md)
  /
  [`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md)).

- `threshold`: the numeric cut-point used.

- `method`: the threshold method applied.

- `train`: the stratified training `data.frame`.

- `histogram`: a `ggplot` object, or `NULL`.

## See also

[`plot_km`](https://fracalanca.github.io/Omics2Surv/reference/plot_km.md),
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md)
