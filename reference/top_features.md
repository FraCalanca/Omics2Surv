# Top Features by Variable Importance

Extracts a model's variable importance from `metadata(mae)$analysis`,
ranks it in decreasing order and cuts it to a feature set. Intended for
models that do *not* perform built-in feature selection (e.g.
`"blockForest"`, `"flexynesis"`); the returned features can then be
passed to
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md).
The model node is located automatically, including methods nested under
container entries such as `"joint"` (e.g. `blockForest` at
`analysis$joint$blockForest`).

## Usage

``` r
top_features(mae, model, n = 50, min_importance = NULL)
```

## Arguments

- mae:

  A `MultiAssayExperiment` whose `metadata()$analysis` holds the fitted
  model.

- model:

  Character string naming the model (e.g. `"blockForest"`,
  `"flexynesis"`).

- n:

  Integer. Number of top features to keep. Default is 50. Set to `NULL`
  to keep all features passing `min_importance`.

- min_importance:

  Numeric. If supplied, only features with importance greater than or
  equal to this value are kept (applied before `n`). Default is `NULL`.

## Value

A character vector of feature names, ordered by decreasing importance.
Importance values are attached as the `"importance"` attribute.

## Details

Both storage formats are handled: a named numeric importance vector, and
a `data.frame` (either the `flexynesis` layout with `name` and
`importance` columns, or a two-column feature/importance table).
Duplicated feature names (e.g. a feature reported across several
`flexynesis` layers) are collapsed by their maximum importance.

## See also

[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md)
