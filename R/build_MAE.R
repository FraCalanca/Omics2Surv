#' Build a MultiAssayExperiment Object from Omics Data
#'
#' Constructs a MultiAssayExperiment (MAE) container that integrates multiple omics
#' datasets with associated clinical/phenotypic data. This unified data structure
#' facilitates multi-omics analysis and ensures coordinated sample management across
#' different data types.
#'
#' @param omics_list A named list of SummarizedExperiment objects, where each element
#'   represents a different omics dataset (e.g., gene expression, methylation, proteomics).
#'   All elements must be valid SummarizedExperiment objects. List names will be used
#'   as experiment identifiers in the resulting MAE.
#' @param clinical_data A matrix or data frame of clinical/phenotypic data where:
#'   \itemize{
#'     \item Rows represent different clinical variables (features)
#'     \item Columns represent samples
#'     \item Column names should match sample identifiers across omics datasets
#'   }
#'   This will be transposed internally so that rows represent samples in the final MAE.
#'
#' @return A MultiAssayExperiment object containing:
#'   \itemize{
#'     \item \code{experiments}: The omics datasets from \code{omics_list}
#'     \item \code{colData}: Transposed and type-converted clinical data as a DataFrame,
#'       where rows represent samples and columns represent clinical variables
#'     \item Coordinated sample mapping across all experiments
#'   }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Validates that all elements in \code{omics_list} are SummarizedExperiment objects
#'   \item Transposes \code{clinical_data} so samples are in rows (standard MAE format)
#'   \item Preserves sample identifiers from rownames after transposition
#'   \item Automatically converts column types using \code{type.convert} (e.g., strings to factors,
#'     numeric strings to numbers) with \code{as.is = FALSE}
#'   \item Constructs the MultiAssayExperiment with proper S4 DataFrame structure
#'   \item Prints a confirmation message with the number of omics datasets included
#' }
#'
#' The resulting MAE object enables:
#' \itemize{
#'   \item Unified access to multiple omics layers
#'   \item Automatic subsetting across all experiments
#'   \item Integration with Bioconductor workflows
#'   \item Sample-level metadata management
#' }
#'
#' @section Data Structure Requirements:
#' \itemize{
#'   \item Each SummarizedExperiment in \code{omics_list} should have sample identifiers
#'     as column names (in \code{colData})
#'   \item Clinical data column names should correspond to these sample identifiers
#'   \item Sample matching across experiments is handled automatically by MAE
#' }
#'
#' @seealso
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} for the MAE class structure
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for omics data structure
#'
#' @export
build_MAE= function (omics_list, clinical_data) {

  if (!all(sapply(omics_list, function(x) inherits(x, "SummarizedExperiment")))) {
    stop("Error: All elements in omics_list must be SummarizedExperiment objects.")
  }

  clinical_data = as.data.frame(t(clinical_data))
  clinical_data$SampleID = rownames(clinical_data)

  clinical_data = as.data.frame(lapply(clinical_data, type.convert, as.is = FALSE))
  rownames(clinical_data) = clinical_data$SampleID
  clinical_data$SampleID = NULL

  mae = MultiAssayExperiment::MultiAssayExperiment(
    experiments = omics_list,
    colData = S4Vectors::DataFrame(clinical_data)
  )

  message("MultiAssayExperiment successfully created with ", length(omics_list), " omics datasets.")

  return(mae)

}
