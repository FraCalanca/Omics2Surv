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


#' Retrieve LinkedOmics Data from TCGA
#'
#' Downloads omics or clinical data for a specific cancer type from LinkedOmics
#' via TCGAbiolinks, performs quality checks, and formats the data for analysis.
#' For omics data, features (rows) with any \code{NA} across samples are removed
#' at this stage so that all downstream functions receive complete matrices.
#'
#' @param cancer_type Character string specifying the TCGA cancer type
#'   (e.g. \code{"TCGA-BRCA"}, \code{"TCGA-LGG"}).
#' @param data_type Character string specifying the data type to retrieve
#'   (e.g. \code{"RNASeq (HiSeq, Gene Level)"},
#'   \code{"Methylation (Gene level, HM450K)"},
#'   \code{"Clinical"}). See \code{TCGAbiolinks::getLinkedOmicsData} for
#'   available types.
#'
#' @return A matrix with:
#'   \itemize{
#'     \item Rows: features (genes, CpG sites, clinical variables, etc.)
#'     \item Columns: TCGA sample identifiers
#'     \item For omics data: features with any \code{NA} are removed.
#'     \item For clinical data: samples with \code{NA} in
#'       \code{"overall_survival"} or \code{"status"} are removed.
#'     \item For methylation data: beta-values are shifted by +0.5 before
#'       return.
#'   }
#'
#' @examples
#' \dontrun{
#' LGG_RNA  <- retrieve_linkedomics_data("TCGA-LGG", "RNASeq (HiSeq, Gene Level)")
#' LGG_Clin <- retrieve_linkedomics_data("TCGA-LGG", "Clinical")
#' }
#'
#' @export
retrieve_linkedomics_data <- function(cancer_type, data_type) {

  linkedomics_data <- as.data.frame(
    TCGAbiolinks::getLinkedOmicsData(cancer_type, data_type)
  )

  # ── Set rownames from attrib_name ────────────────────────────────
  if (!"attrib_name" %in% colnames(linkedomics_data)) {
    stop("No 'attrib_name' column found in the downloaded data.")
  }
  rownames(linkedomics_data)  <- linkedomics_data$attrib_name
  linkedomics_data$attrib_name <- NULL

  # ── Warn on duplicated sample IDs ────────────────────────────────
  if (length(unique(colnames(linkedomics_data))) != ncol(linkedomics_data)) {
    warning("Duplicated sample identifiers detected.")
  }

  # ── Clinical: remove samples with NA survival ────────────────────
  if (data_type == "Clinical") {

    missing_rows <- setdiff(c("overall_survival", "status"),
                            rownames(linkedomics_data))
    if (length(missing_rows) > 0L) {
      stop("Missing required rows: ",
           paste(missing_rows, collapse = ", "))
    }

    na_samples <- is.na(linkedomics_data["overall_survival", ]) |
      is.na(linkedomics_data["status", ])
    if (any(na_samples)) {
      linkedomics_data <- linkedomics_data[, !na_samples, drop = FALSE]
      message(sum(na_samples), " sample(s) with NA in 'overall_survival' or ",
              "'status' removed.")
    }

  } else {

    # ── Omics: remove features (rows) with any NA ─────────────────
    na_features <- rowSums(is.na(linkedomics_data)) > 0L
    if (any(na_features)) {
      linkedomics_data <- linkedomics_data[!na_features, , drop = FALSE]
      message(sum(na_features), " feature(s) with NA removed. ",
              nrow(linkedomics_data), " features remaining.")
    }

    # ── Methylation: shift beta-values ────────────────────────────
    methyl_types <- c("Methylation (Gene level, HM27)",
                      "Methylation (Gene level, HM450K)")
    if (data_type %in% methyl_types) {
      linkedomics_data <- linkedomics_data + 0.5
    }
  }

  as.matrix(linkedomics_data)
}


#' Prepare Custom Data in LinkedOmics-like Format
#'
#' Validates and cleans a user-supplied omics or clinical matrix, applying the
#' same quality checks and NA handling used by
#' \code{\link{retrieve_linkedomics_data}}. Use this function when data come
#' from a source other than LinkedOmics but must conform to the same format
#' expected by the rest of the pipeline.
#'
#' @param data A data frame or matrix in LinkedOmics-like format: features as
#'   rows, samples as columns, with row names set to feature identifiers.
#' @param data_type Character string specifying the data type:
#'   \itemize{
#'     \item \code{"Omics"}: features with any \code{NA} are removed.
#'     \item \code{"Clinical"}: samples with \code{NA} in
#'       \code{"overall_survival"} or \code{"status"} are removed.
#'   }
#'
#' @return A matrix with:
#'   \itemize{
#'     \item For \code{"Omics"}: features (rows) with any \code{NA} removed.
#'     \item For \code{"Clinical"}: samples (columns) with \code{NA} in
#'       survival variables removed.
#'   }
#'
#' @seealso \code{\link{retrieve_linkedomics_data}}
#'
#' @export
prepare_custom_data <- function(data,
                                data_type = c("Omics", "Clinical")) {

  data_type <- match.arg(data_type)
  data      <- as.matrix(data)

  # ── Validate row names ───────────────────────────────────────────
  if (is.null(rownames(data))) {
    stop("Feature identifiers must be provided as row names.")
  }

  # ── Warn on duplicated sample IDs ────────────────────────────────
  if (anyDuplicated(colnames(data))) {
    warning("Duplicated sample identifiers detected.")
  }

  # ── Type-specific cleaning ────────────────────────────────────────
  if (data_type == "Clinical") {

    required_rows <- c("overall_survival", "status")
    missing_rows  <- setdiff(required_rows, rownames(data))
    if (length(missing_rows) > 0L) {
      stop("Missing required clinical rows: ",
           paste(missing_rows, collapse = ", "))
    }

    na_samples <- is.na(data["overall_survival", ]) |
      is.na(data["status", ])
    if (any(na_samples)) {
      data <- data[, !na_samples, drop = FALSE]
      message(sum(na_samples), " sample(s) with NA in 'overall_survival' or ",
              "'status' removed.")
    }

  } else {

    # ── Omics: remove features (rows) with any NA ─────────────────
    na_features <- rowSums(is.na(data)) > 0L
    if (any(na_features)) {
      data <- data[!na_features, , drop = FALSE]
      message(sum(na_features), " feature(s) with NA removed. ",
              nrow(data), " features remaining.")
    }
  }

  data
}


#' Split Multi-Omics Data into Training and Test Sets
#'
#' Identifies the common sample set across all omics layers, performs a single
#' stratified random split on that common set, and returns the raw matrices
#' partitioned into training and test subsets. This ensures every patient lands
#' in the same split across all omics layers, which is a prerequisite for
#' multi-omics survival analysis.
#'
#' @param omics_list Named list of raw omics matrices (features as rows, samples
#'   as columns). Matrices do not need to share the same samples beforehand —
#'   the intersection is computed internally.
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns). Must contain \code{"overall_survival"} and
#'   \code{"status"} rows, and must cover all common samples.
#' @param train_fraction Numeric value in (0, 1) specifying the proportion of
#'   common samples to assign to training. Default is 0.7.
#' @param seed Integer seed for reproducibility. Default is 4820.
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{train_ids}: Character vector of training sample IDs.
#'     \item \code{test_ids}: Character vector of test sample IDs.
#'     \item \code{train_list}: Named list of raw matrices subset to training
#'       samples (same names as \code{omics_list}).
#'     \item \code{test_list}: Named list of raw matrices subset to test samples.
#'     \item \code{train_clinical}: Clinical data matrix subset to training
#'       samples (variables as rows, samples as columns).
#'     \item \code{test_clinical}: Clinical data matrix subset to test samples.
#'   }
#'
#' @details
#' The split is performed once on the intersection of all sample sets and then
#' applied identically to every omic layer and to \code{clinical_data}. This
#' guarantees that:
#' \itemize{
#'   \item No patient appears in training for one omic and test for another.
#'   \item \code{build_MAE} on the resulting train/test lists will not drop
#'     any additional samples when computing the common set.
#' }
#'
#' Samples present in some but not all omics are silently excluded from both
#' splits. A message reports the common sample count and how many samples were
#' excluded per layer.
#'
#' @seealso \code{\link{prepare_RNA}}, \code{\link{prepare_test_SE}},
#'   \code{\link{build_MAE}}
#'
#' @export
split_train_test <- function(
    omics_list,
    clinical_data,
    train_fraction = 0.7,
    seed           = 4820
) {

  # ── Validate inputs ──────────────────────────────────────────────
  if (!is.list(omics_list) || is.null(names(omics_list))) {
    stop("'omics_list' must be a named list of matrices.")
  }
  if (train_fraction <= 0 || train_fraction >= 1) {
    stop("'train_fraction' must be strictly between 0 and 1.")
  }

  required_rows <- c("overall_survival", "status")
  missing_rows  <- setdiff(required_rows, rownames(clinical_data))
  if (length(missing_rows) > 0L) {
    stop("'clinical_data' missing required rows: ",
         paste(missing_rows, collapse = ", "))
  }

  # ── Intersect samples across all omics AND clinical_data ─────────
  common_samples <- Reduce(
    intersect,
    c(lapply(omics_list, colnames), list(colnames(clinical_data)))
  )

  if (length(common_samples) == 0L) {
    stop("No common samples found across all omics layers and clinical_data.")
  }

  n_before <- sapply(omics_list, ncol)
  dropped  <- n_before - length(common_samples)

  message("Common samples across all omics and clinical data: ",
          length(common_samples))
  for (nm in names(dropped)[dropped > 0L]) {
    message("  '", nm, "': ", dropped[nm],
            " sample(s) excluded (not present in all layers or clinical data).")
  }

  # ── Single split on common samples ───────────────────────────────
  set.seed(seed)
  n_total <- length(common_samples)
  n_train <- round(n_total * train_fraction)

  train_idx  <- sample(seq_len(n_total), size = n_train, replace = FALSE)
  test_idx   <- setdiff(seq_len(n_total), train_idx)

  train_ids  <- common_samples[train_idx]
  test_ids   <- common_samples[test_idx]

  message("Training set: ", length(train_ids), " samples.")
  message("Test set:     ", length(test_ids),  " samples.")

  # ── Subset each omic ─────────────────────────────────────────────
  train_list <- lapply(omics_list, function(mat) mat[, train_ids, drop = FALSE])
  test_list  <- lapply(omics_list, function(mat) mat[, test_ids,  drop = FALSE])

  # ── Subset clinical data ─────────────────────────────────────────
  train_clinical <- clinical_data[, train_ids, drop = FALSE]
  test_clinical  <- clinical_data[, test_ids,  drop = FALSE]

  list(
    train_ids      = train_ids,
    test_ids       = test_ids,
    train_list     = train_list,
    test_list      = test_list,
    train_clinical = train_clinical,
    test_clinical  = test_clinical
  )
}
