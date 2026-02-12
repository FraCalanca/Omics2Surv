#' Prepare Design Matrix from MultiAssayExperiment for Modeling
#'
#' Extracts omics data from a MultiAssayExperiment object and optionally combines it
#' with clinical variables to create a design matrix suitable for survival modeling.
#' Handles sample ordering, NA removal, and factor variable conversion.
#'
#' @param mae A MultiAssayExperiment object containing omics experiments and clinical
#'   data in \code{colData}.
#' @param assay_name Character string specifying the name of the assay (omics layer)
#'   to extract from the MAE. Must match an experiment name in \code{names(mae)}.
#' @param clinical_list Character vector of clinical variable names to include in the
#'   design matrix. Variables are extracted from \code{colData(mae)}. If NULL, only
#'   omics features are included. Default is NULL.
#' @param remove_na Logical indicating whether to remove samples (rows) with any
#'   missing values. If TRUE, samples with NA in any variable are excluded. If FALSE,
#'   NAs are retained. Default is TRUE.
#'
#' @return A numeric matrix where:
#'   \itemize{
#'     \item Rows represent samples (patients)
#'     \item Columns represent features:
#'       \itemize{
#'         \item Clinical variables (if \code{clinical_list} is specified), followed by
#'         \item Omics features (genes, proteins, metabolites, etc.)
#'       }
#'     \item Row names are sample identifiers
#'     \item Column names are feature names
#'     \item All values are numeric (factors are converted to numeric codes)
#'   }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item \strong{Extract omics data}: Retrieves the "counts" assay from the specified
#'     experiment in the MAE
#'   \item \strong{Determine sample order}: Uses column names from omics data as the
#'     reference sample ordering
#'   \item \strong{Extract clinical data}: Retrieves clinical variables from MAE colData,
#'     ordered to match omics samples
#'   \item \strong{Combine data} (if clinical_list is specified):
#'     \itemize{
#'       \item Subsets clinical variables to requested list
#'       \item Converts factors to numeric codes
#'       \item Binds clinical variables (columns) with transposed omics data
#'     }
#'   \item \strong{Handle missing data} (if remove_na = TRUE):
#'     \itemize{
#'       \item Identifies samples with any NA values
#'       \item Removes these samples from the design matrix
#'     }
#' }
#'
#' @seealso
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} for the MAE object class
#' \code{\link[MultiAssayExperiment]{colData}} for accessing clinical data
#' \code{\link{prepare_surv_data}} for preparing survival outcome data
#'
#' @export
prepare_design_matrix <- function(
    mae,
    assay_name,
    clinical_list = NULL,
    remove_na = TRUE
) {

  omic_data <- MultiAssayExperiment::assays(mae[[assay_name]])$counts
  sample_order <- colnames(omic_data)

  clinical_data <- as.data.frame(MultiAssayExperiment::colData(mae))[sample_order, , drop = FALSE]

  if (!is.null(clinical_list)) {
    clinical_vars <- clinical_data[, clinical_list, drop = FALSE]
    clinical_vars[] <- lapply(clinical_vars, function(x) {
      if (is.factor(x)) as.numeric(x) else x
    })
    design_matrix <- cbind(clinical_vars, t(omic_data))
  } else {
    design_matrix <- t(omic_data)
  }

  if (remove_na) {
    na_row <- which(rowSums(is.na(design_matrix)) > 0)
    if (length(na_row) > 0) {
      design_matrix <- design_matrix[-na_row, , drop = FALSE]
    }
  }

  design_matrix
}
