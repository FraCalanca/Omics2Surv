# Supervised Feature Screening via Marginal Cox Model

Implements marginal screening for high-dimensional omics features based
on their association with survival outcomes. Features are selected by
one of four criteria or a fixed top-N.

## Usage

``` r
supervised_screening(
  data,
  clinical_data,
  rank = c("pval", "absmg", "padj", "mgpval"),
  p_threshold = 0.05,
  top_n_features = NULL
)
```

## Arguments

- data:

  Numeric matrix of omics data (features as rows, samples as columns).
  Row names must be set.

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns). Must contain rows `"overall_survival"` and `"status"`.

- rank:

  Character string specifying the ranking criterion:

  - `"pval"`: p-value of the marginal Cox test. Supports `p_threshold`
    and `top_n_features`.

  - `"padj"`: Benjamini-Hochberg adjusted p-value. Supports
    `p_threshold` and `top_n_features`.

  - `"absmg"`: absolute marginal Cox coefficient \\\|\hat\beta_j\|\\.
    **Requires `top_n_features`**.

  - `"mgpval"`: combined score \\-\log\_{10}(p) \times
    \|\hat\beta_j\|\\. **Requires `top_n_features`**.

  Default is `"pval"`.

- p_threshold:

  Numeric threshold used only when `rank` is `"pval"` (p-value \\\le\\
  threshold) or `"padj"` (BH-adjusted p-value \\\le\\ threshold), and
  `top_n_features` is `NULL`. Default is `0.05`.

- top_n_features:

  Integer. Number of top-ranked features to retain. **Mandatory** for
  `"absmg"` and `"mgpval"`, which have no natural threshold. Optional
  for `"pval"` and `"padj"` (overrides `p_threshold` when set). Default
  is `NULL`.

## Value

A named list with:

- `data`: Filtered matrix (selected features x samples), with *original*
  (non-standardised) values.

- `rank_table`: Data frame with columns `feature`, `beta`, `p_value`,
  `p_adj`, `score` (\\-\log\_{10}(p) \times \|\beta\|\\), `rank_value`
  (the metric used for sorting), and `rank_position`, for *all*
  features.

- `selected`: Character vector of selected feature names.

- `rank`: The `rank` argument used.

## Details

Designed to be reused across omics `prepare_*` functions (mRNA, miRNA,
CNV, Methylation).
