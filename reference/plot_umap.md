# UMAP of Patients Coloured by Risk Group

Projects patients into two dimensions with
[`umap`](https://jlmelville.github.io/uwot/reference/umap.html) and
colours them by the risk group from
[`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md).
The projection is restricted to the features of the given `model`:
[`selected_features()`](https://fracalanca.github.io/Omics2Surv/reference/selected_features.md)
for models with built-in selection (Cox penalised models) and
[`top_features()`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)
for importance-based models (`"blockForest"`, `"flexynesis"`); the
choice is made from `model`. Assay samples are matched to stratified
patients by barcode.

## Usage

``` r
plot_umap(
  mae,
  x,
  model,
  modality = NULL,
  top_feat = NULL,
  mae_train = NULL,
  features = NULL,
  experiment = NULL,
  assay = 1,
  colors = .risk_colors(),
  rho = NULL,
  scale = TRUE,
  n_neighbors = 15,
  seed = 42,
  ...
)
```

## Arguments

- mae:

  A `MultiAssayExperiment` providing the feature matrix to project
  (typically the stratified test object).

- x:

  Either the `list` returned by
  [`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md)
  or its `strata` `data.frame`.

- model:

  Character string naming the model whose features drive the projection
  (e.g. `"cox_lasso"`, `"blockForest"`, `"flexynesis"`).

- modality:

  Character string naming the omics modality. It locates a penalised
  model's selected features (ignored for `"blockForest"` /
  `"flexynesis"`) and, unless `experiment` is given, also selects the
  assay to project on. Default is `NULL`: the single per-run integration
  tag (`"early_<omics>"`) is then resolved automatically; pass the tag
  (or an omic name) explicitly to disambiguate.

- top_feat:

  Integer number of top features to keep. Used *only* when `model` is
  `"blockForest"` or `"flexynesis"`; ignored otherwise. If `NULL`
  (default),
  [`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)'s
  own default is used.

- mae_train:

  A `MultiAssayExperiment` from which to read the model's features.
  Defaults to `mae`; set this when feature information lives in the
  training object.

- features:

  Optional character vector of feature (row) names. If supplied, it
  overrides the automatic retrieval from `model`.

- experiment:

  Name or index of the experiment to project on. If `NULL` (default),
  the assay matching `modality` is used, or the single `"early_*"`
  integration assay when `modality` is `NULL` (joint models).

- assay:

  Assay index or name within `experiment`. Default is 1.

- colors:

  Named character vector mapping risk groups to colours. Default is
  `c("High risk" = "red", "Low risk" = "deepskyblue")`.

- rho:

  Numeric. For `model = "AFTCoop"` only: selects the `beta_rho_<rho>`
  feature set. If `NULL` (default) and only one rho was trained, it is
  resolved automatically.

- scale:

  Logical. If `TRUE`, features are standardised and zero-variance
  features are dropped. Default is `TRUE`.

- n_neighbors:

  Integer UMAP neighbourhood size. Default is 15.

- seed:

  Integer random seed for reproducibility. Default is 42.

- ...:

  Additional arguments passed to
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

A `ggplot` object.

## Details

Feature information (selected features / variable importance) is a
training artifact, so it is read from `mae_train`. The projection itself
is drawn on `mae` (typically the stratified test set).

## See also

[`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md),
[`plot_km`](https://fracalanca.github.io/Omics2Surv/reference/plot_km.md),
[`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md),
[`selected_features`](https://fracalanca.github.io/Omics2Surv/reference/selected_features.md)
