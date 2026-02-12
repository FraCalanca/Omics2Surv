# Filter Omics Data to Retain Only Common Samples Across All Layers

Identifies and retains only samples (patients) that are present across
all omics datasets. This ensures complete multi-omics profiles for
downstream integrative analysis by removing samples with missing data in
any omics layer.

## Usage

``` r
filter_common_patients(omics_list)
```

## Arguments

- omics_list:

  A list of omics data objects (typically SummarizedExperiment objects
  or matrices) where:

  - Each element represents a different omics layer (e.g.,
    transcriptomics, proteomics, methylation)

  - Samples are represented as columns

  - Column names contain sample identifiers that should be consistent
    across layers

## Value

A list of the same length and structure as `omics_list`, where each
element has been subset to include only the common samples. The order of
samples (columns) is consistent across all returned omics layers. Sample
order follows the intersection order determined by
`Reduce(intersect, ...)`.

## Details

The function performs the following operations:

1.  Extracts column names (sample IDs) from each omics layer

2.  Identifies samples present in ALL layers using set intersection

3.  Subsets each omics object to retain only common samples

4.  Preserves the original data structure and class of each omics layer

This filtering step is essential for multi-omics integration because:

- Many integration methods require complete data across all modalities

- It prevents errors from sample mismatches during analysis

- It ensures fair comparison across omics layers

- It identifies the cohort with comprehensive multi-omics profiling

The function will stop with an error if no samples are common across all
layers, which may indicate:

- Inconsistent sample naming conventions

- Non-overlapping patient cohorts
