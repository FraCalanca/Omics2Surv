# Dotplot of an Enrichment Result

Convenience wrapper around
[`dotplot`](https://rdrr.io/pkg/enrichplot/man/dotplot.html) for one
database of an
[`enrich_features`](https://fracalanca.github.io/Omics2Surv/reference/enrich_features.md)
result.

## Usage

``` r
plot_enrichment(
  enrichment,
  database = c("GO", "KEGG", "Reactome"),
  show_category = 15,
  ...
)
```

## Arguments

- enrichment:

  The list returned by
  [`enrich_features`](https://fracalanca.github.io/Omics2Surv/reference/enrich_features.md),
  or a single `enrichResult` object.

- database:

  Which result to plot when `enrichment` is a list: `"GO"`, `"KEGG"` or
  `"Reactome"`. Default is `"GO"`.

- show_category:

  Integer number of top terms to display. Default is 15.

- ...:

  Additional arguments passed to
  [`dotplot`](https://rdrr.io/pkg/enrichplot/man/dotplot.html).

## Value

A `ggplot` object.

## See also

[`enrich_features`](https://fracalanca.github.io/Omics2Surv/reference/enrich_features.md)
