# Kaplan-Meier Plot by Risk Group

Draws Kaplan-Meier survival curves for the `"High risk"` and
`"Low risk"` groups produced by
[`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md),
using [`ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html).

## Usage

``` r
plot_km(
  x,
  colors = .risk_colors(),
  pval = TRUE,
  conf_int = TRUE,
  risk_table = TRUE,
  time_unit = "Time",
  ...
)
```

## Arguments

- x:

  Either the `list` returned by
  [`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md)
  or its `strata` `data.frame`.

- colors:

  Named character vector mapping risk groups to colours. Default is
  `c("High risk" = "red", "Low risk" = "deepskyblue")`.

- pval:

  Logical, passed to
  [`ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html).
  Default is `TRUE`.

- conf_int:

  Logical, passed to
  [`ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html) as
  `conf.int` to draw the confidence bands. Default is `TRUE`.

- risk_table:

  Logical, passed to
  [`ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html) as
  `risk.table`. Default is `TRUE`.

- time_unit:

  Character string used as the x-axis label. Default is `"Time"`.

- ...:

  Additional arguments passed to
  [`ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html).

## Value

A `ggsurvplot` object.

## See also

[`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md),
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md)
