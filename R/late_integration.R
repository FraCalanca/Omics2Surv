#' Late Integration of Multi-Omics Survival Predictions
#'
#' Combines risk predictions from multiple independently-trained omics models into
#' a single consensus prediction. Supports averaging of normalized ranks or Robust
#' Rank Aggregation (RRA) to create an ensemble prediction that leverages information
#' from all omics layers.
#'
#' @param res_list A named list where each element contains results from a single
#'   omics layer model. Each element must include:
#'   \itemize{
#'     \item \code{pred}: Numeric vector of risk predictions for samples
#'   }
#'   List names typically represent omics layers (e.g., "transcriptomics", "proteomics").
#' @param y A data frame containing survival outcomes with columns:
#'   \itemize{
#'     \item \code{overall_survival}: Numeric survival time
#'     \item \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#'   Row names should correspond to patient/sample identifiers.
#' @param model Character string specifying the underlying prediction model type.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"cox_lasso"}: Cox Lasso regression
#'     \item \code{"cox_adaptive_lasso"}: Adaptive Cox Lasso
#'     \item \code{"cox_elastic_net"}: Cox elastic net
#'     \item \code{"aft_pen"}: Penalized Accelerated Failure Time model
#'   }
#'
#' @param method Character string specifying the integration method. Must be one of:
#'   \itemize{
#'     \item \code{"average"}: Average of normalized ranks across omics layers
#'     \item \code{"RRA"}: Robust Rank Aggregation using the RobustRankAggreg package
#'   }
#'
#' @return A list containing:
#'   \item{pred}{Numeric vector of consensus risk ranks for all patients. Lower ranks
#'     indicate lower risk (longer predicted survival). Length equals number of patients.}
#'   \item{merged_risks}{Data frame containing:
#'     \itemize{
#'       \item \code{patient_id}: Sample identifiers
#'       \item \code{risk_omics1}, \code{risk_omics2}, ...: Original risk scores from each omics layer
#'       \item \code{norm_omics1}, \code{norm_omics2}, ...: Z-score normalized risk scores
#'       \item \code{rank_omics1}, \code{rank_omics2}, ...: Ranks of normalized scores within each layer
#'       \item \code{consensus_rank}: Final integrated consensus rank
#'     }}
#'   \item{cindex}{Numeric concordance index (C-index) evaluating the consensus predictions
#'     against the observed survival outcomes.}
#'
#' @details
#' The late integration workflow consists of the following steps:
#' \enumerate{
#'   \item \strong{Extract}: Collect risk scores from each omics layer model
#'   \item \strong{Merge}: Combine risk scores into a single data frame by patient ID
#'   \item \strong{Normalize}: Z-score standardize risk scores within each omics layer
#'     (mean = 0, SD = 1) to make them comparable across layers
#'   \item \strong{Rank}: Convert normalized scores to ranks within each layer using
#'     average method for ties
#'   \item \strong{Integrate}: Combine ranks across layers:
#'     \itemize{
#'       \item \code{average}: Simple arithmetic mean of ranks
#'       \item \code{RRA}: Statistical aggregation accounting for rank distribution
#'     }
#'   \item \strong{Evaluate}: Calculate C-index of consensus predictions
#' }
#'
#' @seealso
#' \code{\link{L_pen_regression_models}} for the wrapper function using late integration
#' \code{\link[RobustRankAggreg]{aggregateRanks}} for RRA implementation
#' \code{\link[survival]{concordance}} for C-index calculation
#'
#' @export
late_integration <- function(res_list, y, model=c("cox_lasso","cox_adaptive_lasso","cox_elastic_net",
                                                           "aft_pen"),
                             method = c("average","RRA")) {

  method <- match.arg(method)
  patient_ids <- rownames(y)

  # 1. risk score
  risk_dfs <- lapply(seq_along(res_list), function(i) {
    data.frame(
      patient_id = patient_ids,
      risk_score = as.numeric(res_list[[i]]$pred)
    )
  })
  names(risk_dfs) <- paste0("omics", seq_along(res_list))

  # 2. merge
  merged_risks <- Reduce(function(d1, d2) merge(d1, d2, by = "patient_id"),
                         lapply(names(risk_dfs), function(nm) {
                           df <- risk_dfs[[nm]]
                           colnames(df)[2] <- paste0("risk_", nm)
                           df
                         }))

  # 3. norm
  for (nm in names(risk_dfs)) {
    col_risk <- paste0("risk_", nm)
    col_norm <- paste0("norm_", nm)
    merged_risks[[col_norm]] <- as.numeric(scale(merged_risks[[col_risk]]))
  }

  # 3b. RANKING
  rank_cols <- c()
  for (nm in names(risk_dfs)) {
    col_norm <- paste0("norm_", nm)
    col_rank <- paste0("rank_", nm)
    merged_risks[[col_rank]] <- rank(merged_risks[[col_norm]], ties.method = "average")
    rank_cols <- c(rank_cols, col_rank)
  }

  # 4. consensus
  if (method == "average") {
    merged_risks$consensus_rank <- rowMeans(merged_risks[, rank_cols])
  } else if (method == "RRA") {
    rank_matrix <- as.matrix(merged_risks[, rank_cols])
    rownames(rank_matrix) <- merged_risks$patient_id

    library(RobustRankAggreg)
    rra_result <- aggregateRanks(rmat = rank_matrix)

    # RRA in consensus rank
    merged_risks$consensus_rank <- match(merged_risks$patient_id, rra_result$Name)
  }

  # 5. C-index
  y_surv <- survival::Surv(y$overall_survival, y$status)

  if(model=="aft_pen"){
    cindex_obj <- survival::concordance(y_surv ~ merged_risks$consensus_rank, reverse = TRUE)
  } else {
    cindex_obj <- survival::concordance(y_surv ~ merged_risks$consensus_rank)
  }

  return(list(
    pred = merged_risks$consensus_rank,
    merged_risks = merged_risks,
    cindex = cindex_obj$concordance
  ))
}
