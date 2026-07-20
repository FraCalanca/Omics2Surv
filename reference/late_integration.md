# Late Integration of Multi-Omics Survival Predictions

Combines risk predictions from multiple independently-trained omics
models into a single consensus prediction. Supports averaging of
normalized ranks or Robust Rank Aggregation (RRA) to create an ensemble
prediction that leverages information from all omics layers.

## Usage

``` r
late_integration(res_list, y, model = NULL, method = c("average", "RRA"))
```

## Arguments

- res_list:

  A named list where each element contains results from a single omics
  layer model. Each element must include:

  - `pred`: Numeric vector of risk predictions for samples

  List names typically represent omics layers (e.g., "transcriptomics",
  "proteomics").

- y:

  A data frame containing survival outcomes with columns:

  - `overall_survival`: Numeric survival time

  - `status`: Numeric event indicator (1 = event, 0 = censored)

  Row names should correspond to patient/sample identifiers.

- model:

  Character string indicating whether the underlying prediction model is
  a Penalized Accelerated Failure Time model (`"aft_pen"`). Default is
  `NULL`, corresponding to the use of Cox models.

- method:

  Character string specifying the integration method. Must be one of:

  - `"average"`: Average of normalized ranks across omics layers

  - `"RRA"`: Robust Rank Aggregation using the RobustRankAggreg package

## Value

A list containing:

- pred:

  Numeric vector of consensus risk ranks for all patients. Lower ranks
  indicate lower risk (longer predicted survival). Length equals number
  of patients.

- merged_risks:

  Data frame containing:

  - `patient_id`: Sample identifiers

  - `risk_omics1`, `risk_omics2`, ...: Original risk scores from each
    omics layer

  - `rank_omics1`, `rank_omics2`, ...: Ranks of normalized scores within
    each layer

  - `consensus_rank`: Final integrated consensus rank

## Details

The late integration workflow consists of the following steps:

1.  **Extract**: Collect risk scores from each omics layer model

2.  **Merge**: Combine risk scores into a single data frame by patient
    ID

3.  **Rank**: Convert normalized scores to ranks within each layer using
    average method for ties

4.  **Integrate**: Combine ranks across layers:

    - `average`: Simple arithmetic mean of ranks

    - `RRA`: Statistical aggregation accounting for rank distribution

## See also

[`aggregateRanks`](https://rdrr.io/pkg/RobustRankAggreg/man/aggregateRanks.html)
for RRA implementation
