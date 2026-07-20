# Exploratory PCA and eigencorplot

Performs PCA on omics data and generates:

1.  2D PCA plot (PC1 vs PC2)

2.  eigencorplot against selected clinical variables

## Usage

``` r
explore_pca(omics_data, clinical_data, metavars)
```

## Arguments

- omics_data:

  Omics data.frame (Linkedomics-like: features × samples)

- clinical_data:

  Clinical data.frame (Linkedomics-like: features × samples)

- metavars:

  Clinical variables for eigencorplot

## Value

list containing:

- pca: PCAtools PCA object

- biplot: 2D PCA plot

- eigencorplot: eigencorplot object
