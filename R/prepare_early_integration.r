#' Early Integration of Multiple Omics Layers
#'
#' Performs early integration (early fusion) of multiple omics datasets by concatenating
#' feature matrices into a single unified matrix. All omics layers are combined before
#' any modeling, creating a comprehensive feature space that includes all measured
#' molecular variables.
#'
#' @param mae A MultiAssayExperiment object containing multiple omics experiments.
#' @param omics_names Character vector specifying the names of omics experiments to
#'   integrate. Must match experiment names in \code{names(mae)}. Typically includes
#'   2-3 layers such as \code{c("transcriptomics", "proteomics", "methylation")}.
#' @param new_layer_name Character string specifying the name for the new integrated
#'   assay that will be added to the MAE. Default is \code{"early_integrated"}.
#'
#' @return A MultiAssayExperiment object containing:
#'   \itemize{
#'     \item All original omics experiments from the input MAE (subset to common samples)
#'     \item A new experiment named \code{new_layer_name} containing the integrated
#'       feature matrix
#'   }
#'   The integrated experiment is a SummarizedExperiment with:
#'   \itemize{
#'     \item Rows: All features from all omics layers (concatenated)
#'     \item Columns: Samples common across all specified omics layers
#'     \item Row names: Prefixed with omics layer name (e.g., "transcriptomics_GENE1")
#'   }
#'
#' @details
#' The early integration workflow consists of:
#' \enumerate{
#'   \item \strong{Sample alignment}: Identifies and retains only samples present in
#'     ALL specified omics layers using \code{intersectColumns()}
#'   \item \strong{Feature concatenation}:
#'     \itemize{
#'       \item Extracts feature matrix from each omics layer
#'       \item Prefixes row names with omics layer identifier to ensure uniqueness
#'       \item Concatenates all matrices by rows (vertical stacking)
#'     }
#'   \item \strong{SummarizedExperiment creation}: Packages the combined matrix into
#'     a SummarizedExperiment object
#'   \item \strong{MAE integration}: Adds the integrated layer to the original MAE
#'     while preserving all original experiments
#' }
#'
#' @export
early_integration <- function(mae, omics_names, new_layer_name = "early_integrated") {

  # 1.
  mae_matched <- intersectColumns(mae[, , omics_names])

  # 2.
  list_of_matrices <- lapply(omics_names, function(name) {
    mat <- assay(mae_matched[[name]])
    rownames(mat) <- paste0(name, "_", rownames(mat))
    return(mat)
  })

  combined_matrix <- do.call(rbind, list_of_matrices)

  # 3. SummarizedExperiment
  combined_se <- SummarizedExperiment(assays = list(counts = combined_matrix))

  # 4.
  final_mae <- c(mae_matched,
                 setNames(list(combined_se), new_layer_name),
                 mapFrom = omics_names[1])

  return(final_mae)
}
