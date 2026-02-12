#' Prepare custom data in LinkedOmics-like format
#'
#' @param data A data.frame or matrix in LinkedOmics-like format
#' @param data_type Character string specifying the data type.
#'        Supported values are "Omics" and "Clinical".
#'
#' @seealso LinkedOmicsLikeData
#'
#' @return A cleaned data.frame ready for downstream analyses.
#'
#' @export
prepare_custom_data = function(data, data_type = c("Omics", "Clinical")) {

  data_type <- match.arg(data_type)

  ## Coerce to data.frame
  data <- as.data.frame(data)

  ## Check rownames
  if (is.null(rownames(data))) {
    stop("Error: Feature identifiers must be provided as row names.")
  }

  ## Check duplicated samples
  if (anyDuplicated(colnames(data))) {
    warning("Duplicated sample names detected.")
  }

  ## Warn on missing values
  if (anyNA(data)) {
    message("Warning: Missing values detected in the data.")
  }

  ## Clinical-specific cleaning
  if (data_type == "Clinical") {

    required_rows <- c("overall_survival", "status")
    missing_rows <- setdiff(required_rows, rownames(data))

    if (length(missing_rows) > 0) {
      stop(
        "Error: Missing required clinical rows: ",
        paste(missing_rows, collapse = ", ")
      )
    }

    ## Remove samples with NA in survival or status
    na_samples <- which(
      is.na(data["overall_survival", ]) |
        is.na(data["status", ])
    )

    if (length(na_samples) > 0) {
      data <- data[, -na_samples, drop = FALSE]
      message(
        "Removed ", length(na_samples),
        " samples with NA in overall_survival or status."
      )
    }
  }

  return(data)
}
