# Feature Heatmap Split by Risk Group

Draws a [`pheatmap`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html) of
a model's features across the stratified patients. Features are selected
exactly as in
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md):
[`selected_features()`](https://fracalanca.github.io/Omics2Surv/reference/selected_features.md)
for penalised models and
[`top_features()`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)
for importance-based models (`"blockForest"`, `"flexynesis"`), chosen
from `model`.

## Usage

``` r
plot_heatmap(
  mae,
  x,
  model,
  modality = NULL,
  top_feat = NULL,
  mae_train = NULL,
  features = NULL,
  rho = NULL,
  experiment = NULL,
  assay = 1,
  colors = .risk_colors(),
  annotations = NULL,
  discrete_max = 5,
  heat_colors = (grDevices::colorRampPalette(c("navy", "white", "firebrick3")))(200),
  cap = 3,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  ...
)
```

## Arguments

- mae:

  A `MultiAssayExperiment` providing the expression matrix (typically
  the stratified test object).

- x:

  Either the `list` returned by
  [`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md)
  or its `strata` `data.frame`.

- model:

  Character string naming the model whose features are shown.

- modality:

  Character string naming the omics modality. Locates a penalised
  model's selected features (ignored for `"blockForest"` /
  `"flexynesis"`) and, unless `experiment` is given, selects the assay
  to draw. Default is `NULL`: the single per-run integration tag
  (`"early_<omics>"`) is resolved automatically; pass the tag (or an
  omic name) explicitly to disambiguate.

- top_feat:

  Integer number of top features. Used *only* for `"blockForest"` /
  `"flexynesis"`. If `NULL` (default),
  [`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)'s
  own default is used.

- mae_train:

  A `MultiAssayExperiment` from which to read the model's features.
  Defaults to `mae`.

- features:

  Optional character vector of feature names, overriding the automatic
  retrieval from `model`.

- rho:

  Numeric. For `model = "AFTCoop"` only: selects the `beta_rho_<rho>`
  feature set. If `NULL` (default) and only one rho was trained, it is
  resolved automatically.

- experiment:

  Name or index of the experiment to draw. If `NULL` (default), the
  assay matching `modality` is used, or the single `"early_*"`
  integration assay when `modality` is `NULL` (joint models).

- assay:

  Assay index or name within `experiment`. Default is 1.

- colors:

  Named character vector mapping risk groups to annotation colours.
  Default is `c("High risk" = "red", "Low risk" = "deepskyblue")`.

- annotations:

  Optional character vector of `colData(mae)` column names to add as
  extra column annotations alongside `Risk` (e.g.
  `c("status", "pathologic_stage")`). Colours for these are assigned by
  `pheatmap`.

- discrete_max:

  Integer. Extra annotation columns are shown as discrete categories
  unless they are numeric with more than `discrete_max` distinct values,
  which stay a continuous gradient. Default is 5.

- heat_colors:

  Character vector of colours for the heatmap body. Default is a
  navy-white-firebrick ramp of 200 steps.

- cap:

  Numeric. Row z-scores are clipped to `[-cap, cap]`. Default is 3.

- cluster_cols:

  Logical. If `TRUE` (default), patients are clustered and their order
  is left to `pheatmap`; if `FALSE`, they are ordered by risk group with
  a gap marking the split.

- cluster_rows:

  Logical, whether to cluster features. Default is `TRUE`.

- clustering_method:

  Agglomeration method for row clustering. Default is `"ward.D2"`.

- show_rownames, show_colnames:

  Logical, whether to print feature / sample labels. Defaults are `TRUE`
  and `FALSE`.

- fontsize_row:

  Font size for feature labels. Default is 6.

- ...:

  Additional arguments passed to
  [`pheatmap`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html).

## Value

The `pheatmap` object (invisibly), drawn as a side effect.

## Details

By default the columns (patients) are clustered, so the risk-group split
is *not* imposed: the `Risk` annotation bar then shows whether the
expression-based clusters line up with the risk groups. Set
`cluster_cols = FALSE` to instead order patients by risk group
(`"High risk"` first) with a column gap marking the split. Rows
(features) are clustered with Ward's method. Values are row z-scored and
winsorised to `[-cap, cap]` so a few outliers do not wash out the
palette. Feature information is read from `mae_train` (a training
artifact); the matrix shown is drawn from `mae`.

## See also

[`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md),
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md),
[`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md),
[`selected_features`](https://fracalanca.github.io/Omics2Surv/reference/selected_features.md)
