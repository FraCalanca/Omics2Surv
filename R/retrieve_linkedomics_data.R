#' Retrieve LinkedOmics Data from TCGA
#'
#' Downloads omics or clinical data for a specific cancer type from LinkedOmics via
#' TCGAbiolinks, performs quality checks, and formats the data for analysis.
#'
#' @param cancer_type Character string specifying the TCGA cancer type abbreviation
#'   (e.g., "BRCA", "LUAD", "COAD").
#' @param data_type Character string specifying the data type to retrieve (e.g.,
#'   "Clinical", "mRNA", "Protein", "miRNA"). See TCGAbiolinks documentation for
#'   available types.
#'
#' @return A data frame with:
#'   \itemize{
#'     \item Rows: Features (genes, proteins, clinical variables, etc.) from the
#'       \code{attrib_name} column
#'     \item Columns: Samples (TCGA patient identifiers)
#'     \item For clinical data: Samples with NA in "overall_survival" or "status"
#'       are removed
#'   }
#'
#' @export
retrieve_linkedomics_data = function(cancer_type, data_type) {

  linkedomics_data = as.data.frame(TCGAbiolinks::getLinkedOmicsData(cancer_type, data_type))

  if (!"attrib_name" %in% colnames(linkedomics_data)) {
    stop("Error: There is no 'attrib_name' column")
  }

  rownames(linkedomics_data) = linkedomics_data$attrib_name
  linkedomics_data$attrib_name = NULL

  if (any(is.na(linkedomics_data))) {
    message("Warning: there are missing data")
  }

  if (length(unique(colnames(linkedomics_data))) != ncol(linkedomics_data)) {
    message("Warning: there are duplicated samples")
  }

  if (data_type == "Clinical") {

    missing_rows = setdiff(c("overall_survival", "status"), rownames(linkedomics_data))
    if (length(missing_rows) > 0) {
      stop(paste("Error: Missing rows:", paste(missing_rows, collapse = ", ")))
    }

    na_samples = which(is.na(linkedomics_data["overall_survival", ]) | is.na(linkedomics_data["status", ]))

    if (length(na_samples) > 0) {

      linkedomics_data = linkedomics_data[, -na_samples]
      message("Message: Removed samples with NA in 'overall_survival' or 'status' rows")
    }
  }

  return(linkedomics_data)
}
