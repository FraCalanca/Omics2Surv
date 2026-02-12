#' Prepare Copy Number Variation Data as SummarizedExperiment
#'
#' Processes copy number variation (CNV) data by filtering constant genes, optionally
#' performing feature selection using supervised (survival-based) or unsupervised
#' (variability-based) methods, and packaging the results into a SummarizedExperiment
#' object for downstream analysis.
#'
#' @param cnv_data A numeric matrix of CNV data where:
#'   \itemize{
#'     \item Rows represent genes/genomic regions
#'     \item Columns represent samples
#'     \item Values typically represent copy number states (e.g., -2, -1, 0, 1, 2)
#'       or continuous copy number estimates
#'   }
#' @param clinical_data A matrix or data frame of clinical data where:
#'   \itemize{
#'     \item Columns represent samples (matching \code{cnv_data} column names)
#'     \item Rows represent clinical variables
#'     \item Must include rows named \code{"overall_survival"} and \code{"status"}
#'       if using supervised feature selection
#'   }
#' @param variability_threshold Numeric value between 0 and 1 specifying the proportion
#'   of most variable genes to retain when \code{method = "unsupervised"}. For example,
#'   0.25 retains the top 25\% most variable genes. Default is 0.25. Ignored if method
#'   is not "unsupervised".
#' @param method Character string specifying the feature selection method. Options are:
#'   \itemize{
#'     \item \code{"supervised"}: Cox regression-based selection using survival outcomes
#'     \item \code{"unsupervised"}: Variability-based selection (no outcome information)
#'     \item \code{NULL}: No additional filtering beyond constant gene removal (default)
#'   }
#'   Default is NULL.
#' @param p_threshold Numeric p-value threshold for supervised feature selection.
#'   Only genes with Cox regression p-values <= \code{p_threshold} are retained.
#'   Required when \code{method = "supervised"}. Ignored otherwise. Default is NULL.
#' @param top_n_supervised Integer specifying the number of top genes to select if
#'   fewer than 100 genes pass the p-value threshold in supervised selection. Acts as
#'   a fallback to ensure a minimum number of features. Default is NULL.
#'
#' @return A SummarizedExperiment object containing:
#'   \itemize{
#'     \item \code{assays}: List with one element named "counts" containing the filtered
#'       CNV matrix
#'     \item \code{metadata}: List with preprocessing information:
#'       \itemize{
#'         \item \code{omics}: "CNV"
#'         \item \code{method}: The feature selection method used
#'       }
#'   }
#'
#' @details
#' The function performs the following processing steps:
#' \enumerate{
#'   \item \strong{Sample alignment}: Identifies and retains only samples present in
#'     both CNV and clinical data
#'   \item \strong{Constant gene removal}: Removes genes with no variation across samples
#'     (all values identical)
#'   \item \strong{Feature selection} (optional):
#'     \itemize{
#'       \item \strong{Supervised}: Fits univariate Cox regression for each gene, retains
#'         genes with p-value <= \code{p_threshold}. Uses parallel processing for speed.
#'         Falls back to top N genes if < 100 pass threshold.
#'       \item \strong{Unsupervised}: Calculates variability as sum of absolute CNV values,
#'         retains top \code{variability_threshold} proportion of genes
#'     }
#'   \item \strong{SummarizedExperiment creation}: Packages filtered data with metadata
#' }
#'
#' @seealso
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for the return object class
#' \code{\link[survival]{coxph}} for Cox regression in supervised selection
#'
#' @export
prepare_CNV_SE = function(cnv_data, clinical_data,
                          variability_threshold = 0.25,
                          method = NULL,
                          p_threshold = NULL,
                          top_n_supervised=NULL) {

  samples_to_keep = intersect(colnames(cnv_data), colnames(clinical_data))
  if (length(samples_to_keep) > 0) {
    cnv_data = cnv_data[, samples_to_keep, drop = FALSE]
    clinical_data = clinical_data[, samples_to_keep, drop = FALSE]
    message("There are: ", length(samples_to_keep), " CNV data associated with clinical data.")
  }

  unique_counts = apply(cnv_data, 1, function(x) length(unique(x)))
  constant_genes = unique_counts == 1
  filtered_cnv_data = cnv_data[!constant_genes, ]
  message("After removing constant genes, there are: ", nrow(filtered_cnv_data), " genes.")

  if (!is.null(method) && method == "supervised") {

    bpparam = if (.Platform$OS.type == "unix") {
     BiocParallel::MulticoreParam(workers = max(1, parallel::detectCores() - 2))
    } else {
      BiocParallel::SnowParam(workers = max(1, parallel::detectCores() - 2))
    }

    message("Using ", BiocParallel::bpworkers(bpparam), " cores for parallel processing.")

    time = as.numeric(clinical_data["overall_survival", ])
    event = as.numeric(clinical_data["status", ])

    survival_analysis = BiocParallel::bplapply(1:nrow(filtered_cnv_data), function(i) {
      gene_expression = as.numeric(filtered_cnv_data[i, ])
      model = survival::coxph(Surv(time, event) ~ gene_expression)
      summary(model)$coefficients[1, 5]
    }, BPPARAM = bpparam)

    p_values = unlist(survival_analysis)
    gene_names = rownames(filtered_cnv_data)
    selected_genes = gene_names[p_values <= p_threshold]

    if (length(selected_genes) < 100) {
      message("Fewer than 100 genes passed the p threshold. Selecting top ", top_n_supervised, " genes instead.")
      ordered_indices = order(p_values, na.last = NA)
      top_indices = head(ordered_indices, top_n_supervised)
      selected_genes = gene_names[top_indices]
    } else {
      message("Genes after p-value filtering (p <= ", p_threshold, "): ", length(selected_genes))
    }

    filtered_cnv_data = filtered_cnv_data[selected_genes, , drop = FALSE]
    filtered_cnv_data_var = filtered_cnv_data

  } else if (!is.null(method) && method == "unsupervised") {

    variability = rowSums(abs(filtered_cnv_data))
    variability_sorted = sort(variability, decreasing = TRUE)
    genes_to_keep = round(length(variability_sorted) * variability_threshold)
    top_genes = names(variability_sorted)[1:genes_to_keep]

    filtered_cnv_data_var = filtered_cnv_data[top_genes, , drop = FALSE]
    message("After variability filtering, there are: ", nrow(filtered_cnv_data_var), " genes.")

  } else {
    message("No additional filtering applied (only constant gene removal).")
    filtered_cnv_data_var = filtered_cnv_data
  }

  se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = filtered_cnv_data_var))

  S4Vectors::metadata(se)$preprocessing <- list(
    omics = "CNV",
    method = method
  )

  return(se)
}
