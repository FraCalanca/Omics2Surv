# Bar Chart of Category Counts for a Column

Tabulates the values of a single column and draws a labelled bar chart
of the counts with ggplot2. Optionally relabels the categories and sets
custom fill colours.

## Usage

``` r
display_counts(data, column, labels = NULL, fill_colors = NULL, title = NULL)
```

## Arguments

- data:

  A data frame containing the column to summarise.

- column:

  Character string naming the column of `data` to count.

- labels:

  Optional named character vector remapping the raw category values
  (names) to display labels (values). Default is `NULL` (raw values are
  used as-is).

- fill_colors:

  Optional vector of fill colours passed to
  [`scale_fill_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html).
  Default is `NULL` (ggplot2 default palette).

- title:

  Optional plot title. Default is `NULL`, which uses
  `paste("Distribution of", column)`.

## Value

A `ggplot` object.
