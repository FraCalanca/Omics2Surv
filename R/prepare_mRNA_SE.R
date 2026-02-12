#' Prepare mRNA Expression Data as SummarizedExperiment
#'
#' Processes mRNA expression data by aligning samples with clinical data, removing
#' genes with missing values or low expression, and optionally performing supervised
#' or unsupervised feature selection.
#'
#' @param RNA_data Numeric matrix of mRNA expression data (genes as rows, samples
#'   as columns).
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns). Must include "overall_survival" and "status" rows for
#'   supervised selection.
#' @param q_threshold Numeric quantile threshold (0-1) for low expression filtering.
#'   Genes with expression at this quantile equal to 0 are removed. Default is 0.25.
#' @param method Character string for feature selection: \code{"supervised"} (Cox
#'   regression-based), \code{"unsupervised"} (high variable genes), or \code{NULL}
#'   (no selection). Default is NULL.
#' @param p_threshold Numeric p-value threshold for supervised selection. Default is 0.1.
#' @param hvg_threshold Numeric percentile (0-100) for highly variable gene selection
#'   in unsupervised method. Default is 25 (top 25% most variable genes).
#'
#' @return A SummarizedExperiment object with:
#'   \itemize{
#'     \item \code{assays$counts}: Filtered mRNA expression matrix
#'     \item \code{metadata$preprocessing}: List with \code{omics = "mRNA"} and
#'       \code{method} used
#'   }
#'
#' @export
prepare_mRNA_SE = function(RNA_data, clinical_data, q_threshold = 0.25, method=NULL, p_threshold = 0.1, hvg_threshold = 25) {

  samples_to_keep = intersect(colnames(RNA_data), colnames(clinical_data))

  if (length(samples_to_keep) > 0) {
    RNA_data = RNA_data[, samples_to_keep, drop = FALSE]
    clinical_data = clinical_data[, samples_to_keep]
    message("There are: ", length(samples_to_keep), " transcriptomic data associated to clinical data")
  }

  message("Removing genes with NA")
  genes_to_keep = complete.cases(RNA_data)
  filtered_RNA_data = RNA_data[genes_to_keep, , drop = FALSE]
  message("Now there are: ", nrow(filtered_RNA_data), " genes.")

  message("Filtering genes with low expression")
  genes_to_remove = unlist(lapply(1:nrow(filtered_RNA_data), function(i) {
    gene_i = as.numeric(filtered_RNA_data[i, ])
    Q_value = quantile(gene_i, probs = q_threshold)
    return(Q_value == 0)
  }))

  filtered_RNA_data = filtered_RNA_data[!genes_to_remove, ]
  message("Now there are: ", nrow(filtered_RNA_data), " genes.")

  if (is.null(method)) {
    message("No additional filtering applied.")
    se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = filtered_RNA_data))
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

    survival_analysis = BiocParallel::bplapply(1:nrow(filtered_RNA_data), function(i) {
      gene_expression = as.numeric(filtered_RNA_data[i, ])
      model = survival::coxph(Surv(time, event) ~ gene_expression)
      summary(model)$coefficients[1, 5]
    }, BPPARAM = bpparam)

    p_values = unlist(survival_analysis)
    gene_names = rownames(filtered_RNA_data)
    selected_genes = gene_names[p_values <= p_threshold]
    message("Genes after p-value filtering (p <= ", p_threshold, "): ", length(selected_genes))

    filtered_RNA_data = filtered_RNA_data[selected_genes, , drop = FALSE]

  } else if (method == "unsupervised") {

    high_variable = M3C::featurefilter(filtered_RNA_data, method = "A", percentile = hvg_threshold)
    filtered_RNA_data = high_variable$filtered_data
    message("Genes after High Variable Genes filtering : ", nrow(filtered_RNA_data))

  } else {

    stop("Error: Invalid method. Use 'supervised', 'unsupervised', or leave it NULL.")
  }

  se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = filtered_RNA_data))

  S4Vectors::metadata(se)$preprocessing <- list(
    omics = "mRNA",
    method = method
  )

  return(se)
}
