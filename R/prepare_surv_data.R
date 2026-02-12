#' Prepare Survival Data from MultiAssayExperiment
#'
#' Extracts survival outcome variables from a MultiAssayExperiment's colData and
#' optionally reorders samples to match a specified order.
#'
#' @param mae A MultiAssayExperiment object containing clinical data in \code{colData}.
#' @param surv_list Character vector of column names in \code{colData(mae)} representing
#'   survival variables (e.g., \code{c("overall_survival", "status")}).
#' @param sample_order Optional character vector of sample identifiers specifying the
#'   desired row order. If NULL, original order is preserved. Default is NULL.
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item Columns: Survival variables specified in \code{surv_list}
#'     \item Rows: Samples, ordered according to \code{sample_order} if provided
#'     \item Row names: Sample identifiers
#'   }
#'
#' @export
prepare_surv_data <- function(mae, surv_list, sample_order = NULL) {

  surv_data <- as.data.frame(MultiAssayExperiment::colData(mae)[, surv_list])

  if (!is.null(sample_order)) {
    surv_data <- surv_data[sample_order, , drop = FALSE]
  }

  surv_data
}
