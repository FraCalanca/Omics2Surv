# Biological Feature Screening via a User-Supplied Gene List

Retains only features present in a user-supplied `gene_list` (e.g.
pathway members, a literature-derived signature, or a custom panel).
Designed to be reused across omics `prepare_*` functions (mRNA, miRNA,
CNV, Methylation).

## Usage

``` r
biological_screening(data, gene_list)
```

## Arguments

- data:

  Numeric matrix of omics data (features as rows, samples as columns).
  Row names must be set to feature identifiers.

- gene_list:

  Character vector of feature names to retain.

## Value

A named list with:

- `data`: Filtered matrix (matched features x samples).

- `selected`: Character vector of retained feature names (intersection
  of `gene_list` and `rownames(data)`).

- `not_found`: Character vector of names in `gene_list` that were absent
  from `data`.

## Details

The function performs a simple intersection between `rownames(data)` and
`gene_list`. Row order in the returned matrix follows the order of
`gene_list`, not the original row order of `data`. A `warning` is issued
for any names in `gene_list` not found in `data`; the function stops
only if the intersection is empty.
