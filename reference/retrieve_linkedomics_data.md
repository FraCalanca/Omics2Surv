# Retrieve LinkedOmics Data from TCGA

Downloads omics or clinical data for a specific cancer type from
LinkedOmics via TCGAbiolinks, performs quality checks, and formats the
data for analysis. For omics data, features (rows) with any `NA` across
samples are removed at this stage so that all downstream functions
receive complete matrices.

## Usage

``` r
retrieve_linkedomics_data(cancer_type, data_type)
```

## Arguments

- cancer_type:

  Character string specifying the TCGA cancer type (e.g. `"TCGA-BRCA"`,
  `"TCGA-LGG"`).

- data_type:

  Character string specifying the data type to retrieve (e.g.
  `"RNASeq (HiSeq, Gene Level)"`, `"Methylation (Gene level, HM450K)"`,
  `"Clinical"`). See
  [`TCGAbiolinks::getLinkedOmicsData`](https://rdrr.io/pkg/TCGAbiolinks/man/getLinkedOmicsData.html)
  for available types.

## Value

A matrix with:

- Rows: features (genes, CpG sites, clinical variables, etc.)

- Columns: TCGA sample identifiers

- For omics data: features with any `NA` are removed.

- For clinical data: samples with `NA` in `"overall_survival"` or
  `"status"` are removed.

- For methylation data: beta-values are shifted by +0.5 before return.

## Examples

``` r
if (FALSE) { # \dontrun{
LGG_RNA  <- retrieve_linkedomics_data("TCGA-LGG", "RNASeq (HiSeq, Gene Level)")
LGG_Clin <- retrieve_linkedomics_data("TCGA-LGG", "Clinical")
} # }
```
