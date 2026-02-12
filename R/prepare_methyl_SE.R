#' Prepare Methylation Data as SummarizedExperiment
#'
#' Processes methylation data by aligning samples with clinical data, removing
#' incomplete cases, and optionally performing supervised or unsupervised feature
#' selection.
#'
#' @param methyl_data Numeric matrix of methylation data (genes/CpG sites as rows,
#'   samples as columns).
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns). Must include "overall_survival" and "status" rows for
#'   supervised selection.
#' @param system_type Character string specifying OS type ("unix" or other). Currently
#'   unused; parallel backend auto-detected from platform.
#' @param method Character string for feature selection: \code{"supervised"} (Cox
#'   regression-based), \code{"unsupervised"} (IQR-based), or \code{NULL} (no selection).
#'   Default is NULL.
#' @param p_threshold Numeric p-value threshold for supervised selection. Default is 0.1.
#' @param iqr_threshold Numeric proportion (0-1) of top IQR genes to retain for
#'   unsupervised selection. Default is 0.25.
#'
#' @return A SummarizedExperiment object with:
#'   \itemize{
#'     \item \code{assays$counts}: Filtered methylation matrix
#'     \item \code{metadata$preprocessing}: List with \code{omics = "methylation"}
#'       and \code{method} used
#'   }
#'
#' @export
prepare_methyl_SE= function(methyl_data, clinical_data, system_type= "unix", method=NULL, p_threshold= 0.1, iqr_threshold = 0.25 ) {

  samples_to_keep = intersect(colnames(methyl_data), colnames(clinical_data))

  if (length(samples_to_keep) > 0) {
    methyl_data = methyl_data[, samples_to_keep, drop = FALSE]
    clinical_data=clinical_data[,samples_to_keep]
    message("There are: ", length(samples_to_keep), " methylation data associated to clinical data")
  }

  genes_to_keep = complete.cases(methyl_data)

  filtered_methyl_data = methyl_data[genes_to_keep, , drop = FALSE]

  message("Now there are: ", nrow(filtered_methyl_data), " genes.")

  if (is.null(method)) {
    message("No additional filtering applied.")
    se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = filtered_methyl_data))
    return(se)
  }

  if (method == "supervised") {

    bpparam = if (.Platform$OS.type == "unix") {
      BiocParallel::MulticoreParam(workers = max(1, parallel::detectCores() - 2))
    } else {
      BiocParallel::SnowParam(workers = max(1, parallel::detectCores() - 2))
    }

    message("Using ", BiocParallel::bpworkers(bpparam), " cores for parallel processing.")

    time = as.numeric(clinical_data["overall_survival", ])
    event = as.numeric(clinical_data["status", ])

    survival_analysis = BiocParallel::bplapply(1:nrow(filtered_methyl_data), function(i) {
      gene_expression = as.numeric(filtered_methyl_data[i, ])
      model = survival::coxph(Surv(time, event) ~ gene_expression)
      summary(model)$coefficients[1, 5]
    }, BPPARAM = bpparam)

    p_values = unlist(survival_analysis)
    gene_names = rownames(filtered_methyl_data)
    selected_genes = gene_names[p_values <= p_threshold]

    message("Genes after adjusted p-value filtering (p <= ", p_threshold, "): ", length(selected_genes))

    filtered_methyl_data = filtered_methyl_data[selected_genes, , drop = FALSE]

  } else if (method == "unsupervised") {
    iqr_values = apply(filtered_methyl_data, 1, function(x) IQR(x))

    iqr_sorted = sort(iqr_values, decreasing = TRUE)
    top_genes = ceiling(length(iqr_sorted) * iqr_threshold)

    genes_to_select = names(iqr_sorted)[1:top_genes]

    message("Genes after Unsupervised filtering: " , length(genes_to_select))

    filtered_methyl_data = filtered_methyl_data[genes_to_select, , drop = FALSE]
  }

  se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = filtered_methyl_data))


  S4Vectors::metadata(se)$preprocessing <- list(
    omics = "methylation",
    method = method
  )

  return(se)
}
