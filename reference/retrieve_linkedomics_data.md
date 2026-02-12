# Retrieve LinkedOmics Data from TCGA

Downloads omics or clinical data for a specific cancer type from
LinkedOmics via TCGAbiolinks, performs quality checks, and formats the
data for analysis.

## Usage

``` r
retrieve_linkedomics_data(cancer_type, data_type)
```

## Arguments

- cancer_type:

  Character string specifying the TCGA cancer type abbreviation (e.g.,
  "BRCA", "LUAD", "COAD").

- data_type:

  Character string specifying the data type to retrieve (e.g.,
  "Clinical", "mRNA", "Protein", "miRNA"). See TCGAbiolinks
  documentation for available types.

## Value

A data frame with:

- Rows: Features (genes, proteins, clinical variables, etc.) from the
  `attrib_name` column

- Columns: Samples (TCGA patient identifiers)

- For clinical data: Samples with NA in "overall_survival" or "status"
  are removed
