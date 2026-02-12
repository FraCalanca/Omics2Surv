# Early Integration of Multiple Omics Layers

Performs early integration (early fusion) of multiple omics datasets by
concatenating feature matrices into a single unified matrix. All omics
layers are combined before any modeling, creating a comprehensive
feature space that includes all measured molecular variables.

## Usage

``` r
early_integration(mae, omics_names, new_layer_name = "early_integrated")
```

## Arguments

- mae:

  A MultiAssayExperiment object containing multiple omics experiments.

- omics_names:

  Character vector specifying the names of omics experiments to
  integrate. Must match experiment names in `names(mae)`. Typically
  includes 2-3 layers such as
  `c("transcriptomics", "proteomics", "methylation")`.

- new_layer_name:

  Character string specifying the name for the new integrated assay that
  will be added to the MAE. Default is `"early_integrated"`.

## Value

A MultiAssayExperiment object containing:

- All original omics experiments from the input MAE (subset to common
  samples)

- A new experiment named `new_layer_name` containing the integrated
  feature matrix

The integrated experiment is a SummarizedExperiment with:

- Rows: All features from all omics layers (concatenated)

- Columns: Samples common across all specified omics layers

- Row names: Prefixed with omics layer name (e.g.,
  "transcriptomics_GENE1")

## Details

The early integration workflow consists of:

1.  **Sample alignment**: Identifies and retains only samples present in
    ALL specified omics layers using `intersectColumns()`

2.  **Feature concatenation**:

    - Extracts feature matrix from each omics layer

    - Prefixes row names with omics layer identifier to ensure
      uniqueness

    - Concatenates all matrices by rows (vertical stacking)

3.  **SummarizedExperiment creation**: Packages the combined matrix into
    a SummarizedExperiment object

4.  **MAE integration**: Adds the integrated layer to the original MAE
    while preserving all original experiments
