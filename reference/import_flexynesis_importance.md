# Import flexynesis feature-importance scores

Reads the `<prefix>.feature_importance.<explainer>.csv` file(s) produced
during training (one per explainer) and stacks them. Columns:
`target_variable`, `target_class`, `target_class_label`, `layer`, `name`
(feature), `importance`, `explainer`.

## Usage

``` r
import_flexynesis_importance(paths)
```

## Arguments

- paths:

  Character vector of feature-importance CSV paths (e.g.
  `run_flexynesis()$feature_importance`).

## Value

A single data frame, or `NULL` if no files exist.
