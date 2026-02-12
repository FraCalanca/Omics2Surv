#' LinkedOmics-like data format
#'
#' @description
#' This package accepts either data retrieved from LinkedOmics
#' or user-provided data, as long as the input follows the same
#' structure as LinkedOmics data after import.
#'
#' @details
#' The expected input is a feature-by-sample matrix.
#'
#' \strong{General requirements:}
#' \itemize{
#'   \item Input must be a \code{data.frame} or \code{matrix}.
#'   \item Rows represent features (genes, proteins, CpGs, or clinical variables).
#'   \item Columns represent samples.
#'   \item Feature identifiers must be stored as row names.
#'   \item Sample identifiers must be stored as column names.
#' }
#'
#' \strong{Omics data:}
#' \itemize{
#'   \item Each row corresponds to a molecular feature.
#'   \item Values should be numeric.
#' }
#'
#' \strong{Clinical data:}
#' \itemize{
#'   \item Rows correspond to clinical variables.
#'   \item The following rows are required:
#'     \itemize{
#'       \item \code{overall_survival}
#'       \item \code{status}
#'     }
#' }
#'
#' @name LinkedOmicsLikeData
#' @docType data
NULL
