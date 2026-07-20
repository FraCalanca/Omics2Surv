# ══════════════════════════════════════════════════════════════════════════════
# Reusable screening engines
# ══════════════════════════════════════════════════════════════════════════════

#' Biological Feature Screening via a User-Supplied Gene List
#'
#' Retains only features present in a user-supplied \code{gene_list}
#' (e.g. pathway members, a literature-derived signature, or a custom panel).
#' Designed to be reused across omics \code{prepare_*} functions
#' (mRNA, miRNA, CNV, Methylation).
#'
#' @param data Numeric matrix of omics data (features as rows, samples as
#'   columns). Row names must be set to feature identifiers.
#' @param gene_list Character vector of feature names to retain.
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{data}: Filtered matrix (matched features x samples).
#'     \item \code{selected}: Character vector of retained feature names
#'       (intersection of \code{gene_list} and \code{rownames(data)}).
#'     \item \code{not_found}: Character vector of names in \code{gene_list}
#'       that were absent from \code{data}.
#'   }
#'
#' @details
#' The function performs a simple intersection between \code{rownames(data)}
#' and \code{gene_list}. Row order in the returned matrix follows the order
#' of \code{gene_list}, not the original row order of \code{data}.
#' A \code{warning} is issued for any names in \code{gene_list} not found in
#' \code{data}; the function stops only if the intersection is empty.
#' @export
biological_screening <- function(data, gene_list) {

  # ── Validate inputs ──────────────────────────────────────────────
  if (is.null(rownames(data))) {
    stop("'data' must have row names (feature identifiers).")
  }

  if (is.null(gene_list) || length(gene_list) == 0L) {
    stop("'gene_list' must be a non-empty character vector.")
  }

  # ── Intersect ────────────────────────────────────────────────────
  matched   <- intersect(gene_list, rownames(data))   # preserves gene_list order
  not_found <- setdiff(gene_list, rownames(data))

  if (length(not_found) > 0L) {
    warning(
      length(not_found), " feature(s) in 'gene_list' not found in 'data': ",
      paste(head(not_found, 10L), collapse = ", "),
      if (length(not_found) > 10L) " ..." else ""
    )
  }

  if (length(matched) == 0L) {
    stop("None of the features in 'gene_list' were found in 'data'.")
  }

  message("Biological screening: retaining ", length(matched), " / ",
          length(gene_list), " requested features.")

  list(
    data      = data[matched, , drop = FALSE],
    selected  = matched,
    not_found = not_found
  )
}


#' Supervised Feature Screening via Marginal Cox Model
#'
#' Implements marginal screening for high-dimensional omics features based on
#' their association with survival outcomes.
#' Features are selected by one of four criteria or a fixed top-N.
#'
#' Designed to be reused across omics \code{prepare_*} functions
#' (mRNA, miRNA, CNV, Methylation).
#'
#' @param data Numeric matrix of omics data (features as rows, samples as
#'   columns). Row names must be set.
#' @param clinical_data Matrix or data frame of clinical data (variables as
#'   rows, samples as columns). Must contain rows \code{"overall_survival"}
#'   and \code{"status"}.
#' @param rank Character string specifying the ranking criterion:
#'   \itemize{
#'     \item \code{"pval"}: p-value of the marginal Cox test. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"padj"}: Benjamini-Hochberg adjusted p-value. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"absmg"}: absolute marginal Cox coefficient
#'       \eqn{|\hat\beta_j|}. \strong{Requires \code{top_n_features}}.
#'     \item \code{"mgpval"}: combined score \eqn{-\log_{10}(p) \times
#'       |\hat\beta_j|}. \strong{Requires \code{top_n_features}}.
#'   }
#'   Default is \code{"pval"}.
#' @param p_threshold Numeric threshold used only when \code{rank} is
#'   \code{"pval"} (p-value \eqn{\le} threshold) or \code{"padj"}
#'   (BH-adjusted p-value \eqn{\le} threshold), and \code{top_n_features}
#'   is \code{NULL}. Default is \code{0.05}.
#' @param top_n_features Integer. Number of top-ranked features to retain.
#'   \strong{Mandatory} for \code{"absmg"} and \code{"mgpval"}, which have
#'   no natural threshold. Optional for \code{"pval"} and \code{"padj"}
#'   (overrides \code{p_threshold} when set). Default is \code{NULL}.
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{data}: Filtered matrix (selected features x samples),
#'       with \emph{original} (non-standardised) values.
#'     \item \code{rank_table}: Data frame with columns \code{feature},
#'       \code{beta}, \code{p_value}, \code{p_adj}, \code{score}
#'       (\eqn{-\log_{10}(p) \times |\beta|}), \code{rank_value} (the metric
#'       used for sorting), and \code{rank_position}, for \emph{all} features.
#'     \item \code{selected}: Character vector of selected feature names.
#'     \item \code{rank}: The \code{rank} argument used.
#'   }
#'
#' @export
supervised_screening <- function(
    data,
    clinical_data,
    rank           = c("pval", "absmg", "padj", "mgpval"),
    p_threshold    = 0.05,
    top_n_features = NULL
) {

  rank <- match.arg(rank)

  # ── Validate inputs ──────────────────────────────────────────────
  if (is.null(rownames(data))) {
    stop("'data' must have row names.")
  }

  required_rows <- c("overall_survival", "status")
  missing_rows  <- setdiff(required_rows, rownames(clinical_data))
  if (length(missing_rows) > 0) {
    stop("'clinical_data' missing rows: ", paste(missing_rows, collapse = ", "))
  }

  # top_n_features is mandatory for methods without a natural threshold
  if (rank %in% c("absmg", "mgpval") && is.null(top_n_features)) {
    stop(
      "'top_n_features' must be specified when supervised_method = '", rank, "'. ",
      "Use 'pval' or 'padj' if you want threshold-based selection."
    )
  }

  time  <- as.numeric(clinical_data["overall_survival", ])
  event <- as.numeric(clinical_data["status", ])

  if (length(time) != ncol(data)) {
    stop("Mismatch between samples in 'data' and 'clinical_data'.")
  }

  # ── Standardisation (Fan et al. 2010) ───────────────────────────
  X        <- t(data)                        # n x p
  center   <- colMeans(X)
  X_c      <- sweep(X, 2L, center, "-")
  euc_norm <- sqrt(colSums(X_c^2))
  euc_norm[euc_norm == 0] <- 1               # guard against zero-variance
  X_std    <- sweep(X_c, 2L, euc_norm, "/")

  y_surv <- survival::Surv(time, event)

  # ── Parallel backend ─────────────────────────────────────────────
  bpparam <- if (.Platform$OS.type == "unix") {
    BiocParallel::MulticoreParam(workers = max(1L, parallel::detectCores() - 2L))
  } else {
    BiocParallel::SnowParam(workers = max(1L, parallel::detectCores() - 2L))
  }
  message("Supervised screening ('", rank, "') using ",
          BiocParallel::bpworkers(bpparam), " cores.")

  # ── Marginal Cox per feature ─────────────────────────────────────
  n_features  <- nrow(data)
  cox_results <- BiocParallel::bplapply(
    seq_len(n_features),
    function(j) {
      fit      <- survival::coxph(y_surv ~ X_std[, j])
      coef_mat <- summary(fit)$coefficients
      list(beta    = coef_mat[1L, 1L],
           p_value = coef_mat[1L, 5L])
    },
    BPPARAM = bpparam
  )

  betas    <- vapply(cox_results, `[[`, numeric(1), "beta")
  p_values <- vapply(cox_results, `[[`, numeric(1), "p_value")
  p_adj    <- p.adjust(p_values, method = "BH")
  scores   <- -log10(p_values + .Machine$double.eps) * abs(betas)

  names(betas) <- names(p_values) <- names(p_adj) <- names(scores) <-
    rownames(data)

  # ── Ranking metric ───────────────────────────────────────────────
  # pval / padj : ascending  (smaller = better)
  # absmg / mgpval: descending (larger  = better)
  rank_value <- switch(rank,
                       pval   = p_values,
                       absmg  = abs(betas),
                       padj = p_adj,
                       mgpval = scores
  )
  decreasing <- !(rank %in% c("pval", "padj"))
  rank_order <- order(rank_value, decreasing = decreasing)

  rank_table <- data.frame(
    feature       = rownames(data)[rank_order],
    beta          = betas   [rank_order],
    p_value       = p_values[rank_order],
    p_adj         = p_adj   [rank_order],
    score         = scores  [rank_order],
    rank_value    = rank_value[rank_order],
    rank_position = seq_len(n_features),
    row.names     = NULL,
    stringsAsFactors = FALSE
  )

  # ── Feature selection ────────────────────────────────────────────
  if (!is.null(top_n_features)) {

    n_keep   <- min(top_n_features, n_features)
    selected <- rank_table$feature[seq_len(n_keep)]
    message("Top ", n_keep, " features selected by '", rank, "'.")

  } else {
    # Only reachable for "pval" and "padj" (others error above)
    selected <- switch(rank,
                       pval   = rank_table$feature[rank_table$p_value <= p_threshold],
                       padj = rank_table$feature[rank_table$p_adj   <= p_threshold]
    )
    message("Features with ", rank, " <= ", p_threshold,
            ": ", length(selected))
  }

  if (length(selected) == 0L) {
    warning("No features passed screening.")
  }

  # ── Return original (non-standardised) data ──────────────────────
  list(
    data       = data[selected, , drop = FALSE],
    rank_table = rank_table,
    selected   = selected,
    rank       = rank
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# RNA expression wrapper
# ══════════════════════════════════════════════════════════════════════════════

#' Prepare RNA Expression Data as SummarizedExperiment
#'
#' Processes RNA expression data (e.g., mRNA or miRNA) by aligning samples with
#' clinical data, removing features with low expression, and
#' optionally performing supervised, unsupervised, or biological variable screening.
#'
#' @param RNA_data Numeric matrix of RNA expression data (features as rows,
#'   samples as columns). Features may represent mRNAs, miRNAs, or other
#'   transcriptomic measurements.
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns). Must include "overall_survival" and "status" rows for
#'   supervised selection.
#' @param q_threshold Numeric quantile threshold (0-1) for low-expression
#'   filtering. Features with expression at this quantile equal to 0 are removed.
#'   Default is 0.25.
#' @param method Character string for variable screening:
#'   \code{"supervised"}, \code{"unsupervised"}, \code{"biological"},
#'   or \code{NULL} (no selection). Default is \code{NULL}.
#' @param supervised_method Character string specifying the supervised screening
#'   criterion, passed to \code{\link{supervised_screening}}:
#'   \itemize{
#'     \item \code{"pval"}: p-value of the marginal Cox test. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"padj"}: Benjamini-Hochberg adjusted p-value. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"absmg"}: absolute marginal Cox coefficient
#'       \eqn{|\hat\beta_j|}. Requires \code{top_n_features}.
#'     \item \code{"mgpval"}: combined score
#'       \eqn{-\log_{10}(p)\times|\hat\beta_j|}. Requires
#'       \code{top_n_features}.
#'   }
#'   Default is \code{"pval"}.
#' @param unsupervised_method Character string specifying the unsupervised
#'   screening method:
#'   \itemize{
#'     \item \code{"A"}: coefficient of variation.
#'     \item \code{"MAD"}: median absolute deviation.
#'     \item \code{"laplacian"}: Laplacian score.
#'   }
#'   Default is \code{"A"}.
#' @param p_threshold Numeric p-value threshold for supervised selection.
#'   Used only when \code{supervised_method} is \code{"pval"} or
#'   \code{"padj"} and \code{top_n_features} is \code{NULL}.
#'   Default is 0.05.
#' @param hvg_threshold Numeric percentile (0-100) for highly variable feature
#'   selection in unsupervised methods. Used only when
#'   \code{top_n_features} is \code{NULL}. Default is 25 (top 25% most
#'   variable features).
#' @param top_n_features Integer. If provided, overrides threshold-based
#'   filtering. Mandatory when \code{supervised_method} is \code{"absmg"}
#'   or \code{"mgpval"}. Default is \code{NULL}.
#' @param gene_list Character vector of feature names to retain.
#'   Used only when \code{method = "biological"}. Features may correspond
#'   to genes, miRNAs, or other transcriptomic entities.
#'   Default is \code{NULL}.
#' @param n_neighbors Integer. Number of neighbors for Laplacian score
#'   computation (\code{unsupervised_method = "laplacian"}). Recommended
#'   to tune based on sample size, approximately \eqn{\sqrt{n}}.
#'   Default is \code{NULL}.
#'
#' @return A \code{SummarizedExperiment} object containing:
#'   \itemize{
#'     \item \code{assays$counts}: Filtered RNA expression matrix.
#'     \item \code{metadata$preprocessing}: A list containing:
#'       \itemize{
#'         \item \code{omics}: Omics type. Always \code{"RNA"} (covers
#'           mRNA, miRNA, and other transcriptomic measurements).
#'         \item \code{method}: Screening method applied.
#'         \item \code{sub_method}: Specific screening criterion used.
#'         \item \code{n_features}: Number of retained features.
#'       }
#'   }
#'
#' @export
prepare_RNA <- function(
    RNA_data,
    clinical_data,
    q_threshold         = 0.25,
    method              = NULL,
    supervised_method   = c("pval", "absmg", "padj", "mgpval"),
    unsupervised_method = c("A", "MAD", "laplacian"),
    p_threshold         = 0.05,
    hvg_threshold       = 25,
    top_n_features      = NULL,
    gene_list           = NULL,
    n_neighbors         = NULL
) {

  # ── Match method arguments ──────────────────────────────────────────────────
  supervised_method   <- match.arg(supervised_method)
  unsupervised_method <- match.arg(unsupervised_method)

  # ── 1. Align samples ────────────────────────────────────────────────────────
  samples_to_keep <- intersect(colnames(RNA_data), colnames(clinical_data))

  if (length(samples_to_keep) == 0) {
    stop("No overlapping samples between RNA_data and clinical_data.")
  }

  RNA_data      <- RNA_data[, samples_to_keep, drop = FALSE]
  clinical_data <- clinical_data[, samples_to_keep, drop = FALSE]
  message("Samples overlapping with clinical data: ", length(samples_to_keep))

  # ── Low-expression filter ────────────────────────────────────────────────
  message("Filtering features with low expression (quantile threshold = ", q_threshold, ").")
  low_expr <- vapply(
    seq_len(nrow(RNA_data)),
    function(i) quantile(as.numeric(RNA_data[i, ]), probs = q_threshold) == 0,
    logical(1)
  )
  RNA_data <- RNA_data[!low_expr, , drop = FALSE]
  message("features remaining: ", nrow(RNA_data))

  # ── No additional filtering ──────────────────────────────────────────────
  if (is.null(method)) {
    message("No additional filtering applied.")
    se <- SummarizedExperiment::SummarizedExperiment(
      assays   = list(counts = RNA_data),
      metadata = list(preprocessing = list(
        omics   = "RNA",
        method  = NULL,
        n_features = nrow(RNA_data)
      ))
    )
    return(se)
  }

  # ── Dispatch to screening method ─────────────────────────────────────────
  RNA_data <- switch(
    method,
    supervised = {
      result <- supervised_screening(
        data           = RNA_data,
        clinical_data  = clinical_data,
        rank           = supervised_method,
        p_threshold    = p_threshold,
        top_n_features = top_n_features
      )
      result$data
    },
    unsupervised = .rna_screen_unsupervised(
      RNA_data, unsupervised_method, hvg_threshold, top_n_features, n_neighbors
    ),
    biological = biological_screening(RNA_data, gene_list)$data,
    stop("Invalid method. Choose 'supervised', 'unsupervised', 'biological', or NULL.")
  )

  # ──Build SummarizedExperiment ───────────────────────────────────────────
  sub_method <- switch(
    method,
    supervised   = supervised_method,
    unsupervised = unsupervised_method,
    biological   = "gene_list"
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays   = list(counts = RNA_data),
    metadata = list(preprocessing = list(
      omics      = "RNA",
      method     = method,
      sub_method = sub_method,
      n_features    = nrow(RNA_data)
    ))
  )

  return(se)
}


#' Unsupervised screening for RNA expression data
#'
#' Internal helper for \code{\link{prepare_RNA}}. Selects highly variable
#' features by coefficient of variation (\code{"A"}), median absolute
#' deviation (\code{"MAD"}), or Laplacian score (\code{"laplacian"}).
#'
#' @param RNA_data Numeric matrix (features x samples).
#' @param unsupervised_method One of \code{"A"}, \code{"MAD"}, \code{"laplacian"}.
#' @param hvg_threshold Percentile (0-100) retained when \code{top_n_features}
#'   is \code{NULL}.
#' @param top_n_features Integer or \code{NULL}. Overrides \code{hvg_threshold}.
#' @param n_neighbors Integer or \code{NULL}. Neighbours for Laplacian score.
#'
#' @return Filtered numeric matrix (features x samples).
#' @keywords internal
#' @noRd
.rna_screen_unsupervised <- function(RNA_data, unsupervised_method,
                                     hvg_threshold, top_n_features, n_neighbors) {

  message("Unsupervised screening method: '", unsupervised_method, "'.")

  if (unsupervised_method %in% c("A", "MAD")) {

    if (!is.null(top_n_features)) {
      pct_to_pass <- (top_n_features / nrow(RNA_data)) * 100
      message("Selecting top ", top_n_features, " features.")
    } else {
      pct_to_pass <- hvg_threshold
      message("Retaining top ", hvg_threshold, "% most variable features.")
    }

    result   <- M3C::featurefilter(RNA_data, method = unsupervised_method,
                                   percentile = pct_to_pass)
    filtered <- result$filtered_data

  } else if (unsupervised_method == "laplacian") {

    if (is.null(n_neighbors)) {
      n_neighbors <- max(2L, round(sqrt(ncol(RNA_data))))
      message("n_neighbors not set; using k = ", n_neighbors,
              " (approx. sqrt(n_samples)).")
    }

    message("Computing Laplacian scores (n_neighbors = ", n_neighbors, ").")

    ls_result  <- Rdimtools::do.lscore(
      X      = t(RNA_data),
      ndim   = nrow(RNA_data),
      nntype = "knn",
      params = list(k = n_neighbors)
    )

    ranked_idx <- ls_result$featidx   # ascending score = best first

    if (!is.null(top_n_features)) {
      n_keep   <- min(top_n_features, length(ranked_idx))
      keep_idx <- ranked_idx[seq_len(n_keep)]
      message("Selecting top ", n_keep, " features by Laplacian score.")
    } else {
      n_keep   <- ceiling((hvg_threshold / 100) * nrow(RNA_data))
      keep_idx <- ranked_idx[seq_len(min(n_keep, length(ranked_idx)))]
      message("Retaining top ", hvg_threshold, "% features by Laplacian score.")
    }

    filtered <- RNA_data[keep_idx, , drop = FALSE]
  }

  message("features after unsupervised filtering: ", nrow(filtered))
  filtered
}


# ══════════════════════════════════════════════════════════════════════════════
# Methylation wrapper
# ══════════════════════════════════════════════════════════════════════════════

#' Prepare Methylation Data as SummarizedExperiment
#'
#' Processes methylation data by aligning samples with clinical data, and optionally performing supervised, unsupervised, or
#' biological feature screening.
#'
#' @param methyl_data Numeric matrix or data.frame of methylation data (features as rows,
#'   samples as columns).
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns). Must include \code{"overall_survival"} and
#'   \code{"status"} rows for supervised selection.
#' @param method Character string for variable screening: \code{"supervised"},
#'   \code{"unsupervised"}, \code{"biological"}, or \code{NULL} (no selection).
#'   Default is \code{NULL}.
#' @param supervised_method Character string specifying the supervised screening
#'   criterion, passed to \code{\link{supervised_screening}}:
#'   \itemize{
#'     \item \code{"pval"}: p-value of the marginal Cox test. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"padj"}: Benjamini-Hochberg adjusted p-value. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"absmg"}: absolute marginal Cox coefficient
#'       \eqn{|\hat\beta_j|}. \strong{Requires \code{top_n_features}}.
#'     \item \code{"mgpval"}: combined score \eqn{-\log_{10}(p) \times
#'       |\hat\beta_j|}. \strong{Requires \code{top_n_features}}.
#'   }
#'   Default is \code{"pval"}.
#' @param unsupervised_method Character string specifying the unsupervised
#'   screening method: \code{"IQR"} (interquartile range), \code{"MAD"}
#'   (median absolute deviation), or \code{"laplacian"} (Laplacian score).
#'   Default is \code{"IQR"}.
#' @param p_threshold Numeric p-value threshold for supervised selection. Used
#'   only when \code{supervised_method} is \code{"pval"} or \code{"padj"}
#'   and \code{top_n_features} is \code{NULL}. Default is 0.05.
#' @param hvg_threshold Numeric percentile (0-100) of top features to retain
#'   in unsupervised methods. Used only when \code{top_n_features} is
#'   \code{NULL}. Default is 25 (top 25%).
#' @param top_n_features Integer. If provided, overrides threshold-based
#'   filtering. \strong{Mandatory} when \code{supervised_method} is
#'   \code{"absmg"} or \code{"mgpval"}. Default is \code{NULL}.
#' @param gene_list Character vector of gene names to retain.
#'   Used only when \code{method = "biological"}. Default is \code{NULL}.
#' @param n_neighbors Integer. Number of neighbors for Laplacian score
#'   computation (\code{unsupervised_method = "laplacian"}). Defaults to
#'   \eqn{\approx \sqrt{n}} if \code{NULL}.
#'
#' @return A \code{SummarizedExperiment} object with:
#'   \itemize{
#'     \item \code{assays$counts}: Filtered methylation matrix.
#'     \item \code{metadata$preprocessing}: List with \code{omics},
#'       \code{method}, \code{sub_method}, and number of features retained.
#'   }
#'
#' @export
prepare_methyl <- function(
    methyl_data,
    clinical_data,
    method              = NULL,
    supervised_method   = c("pval", "absmg", "padj", "mgpval"),
    unsupervised_method = c("IQR", "MAD", "laplacian"),
    p_threshold         = 0.05,
    hvg_threshold       = 25,
    top_n_features      = NULL,
    gene_list           = NULL,
    n_neighbors         = NULL
) {

  # ── Match method arguments ──────────────────────────────────────────────────
  supervised_method   <- match.arg(supervised_method)
  unsupervised_method <- match.arg(unsupervised_method)

  # ── 1. Align samples ────────────────────────────────────────────────────────
  samples_to_keep <- intersect(colnames(methyl_data), colnames(clinical_data))

  if (length(samples_to_keep) == 0) {
    stop("No overlapping samples between methyl_data and clinical_data.")
  }

  methyl_data   <- methyl_data[, samples_to_keep, drop = FALSE]
  clinical_data <- clinical_data[, samples_to_keep, drop = FALSE]
  message("Samples overlapping with clinical data: ", length(samples_to_keep))

  # ──  No additional filtering ──────────────────────────────────────────────
  if (is.null(method)) {
    message("No additional filtering applied.")
    se <- SummarizedExperiment::SummarizedExperiment(
      assays   = list(counts = methyl_data),
      metadata = list(preprocessing = list(
        omics       = "methylation",
        method      = NULL,
        n_features  = nrow(methyl_data)
      ))
    )
    return(se)
  }

  # ── Dispatch to screening method ─────────────────────────────────────────
  methyl_data <- switch(
    method,
    supervised = {
      result <- supervised_screening(
        data           = methyl_data,
        clinical_data  = clinical_data,
        rank           = supervised_method,
        p_threshold    = p_threshold,
        top_n_features = top_n_features
      )
      result$data
    },
    unsupervised = .methyl_screen_unsupervised(
      methyl_data, unsupervised_method, hvg_threshold, top_n_features, n_neighbors
    ),
    biological = biological_screening(methyl_data, gene_list)$data,
    stop("Invalid method. Choose 'supervised', 'unsupervised', 'biological', or NULL.")
  )

  # ──  Build SummarizedExperiment ───────────────────────────────────────────
  sub_method <- switch(
    method,
    supervised   = supervised_method,
    unsupervised = unsupervised_method,
    biological   = "gene_list"
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays   = list(counts = methyl_data),
    metadata = list(preprocessing = list(
      omics      = "methylation",
      method     = method,
      sub_method = sub_method,
      n_features = nrow(methyl_data)
    ))
  )

  return(se)
}


#' Unsupervised screening for methylation data
#'
#' Internal helper for \code{\link{prepare_methyl}}. Selects highly variable
#' features by interquartile range (\code{"IQR"}), median absolute deviation
#' (\code{"MAD"}), or Laplacian score (\code{"laplacian"}).
#'
#' @param methyl_data Numeric matrix (features x samples).
#' @param unsupervised_method One of \code{"IQR"}, \code{"MAD"}, \code{"laplacian"}.
#' @param hvg_threshold Percentile (0-100) retained when \code{top_n_features}
#'   is \code{NULL}.
#' @param top_n_features Integer or \code{NULL}. Overrides \code{hvg_threshold}.
#' @param n_neighbors Integer or \code{NULL}. Neighbours for Laplacian score.
#'
#' @return Filtered numeric matrix (features x samples).
#' @keywords internal
#' @noRd
.methyl_screen_unsupervised <- function(methyl_data, unsupervised_method,
                                        hvg_threshold, top_n_features, n_neighbors) {

  message("Unsupervised screening method: '", unsupervised_method, "'.")

  if (unsupervised_method %in% c("IQR", "MAD")) {

    metric <- if (unsupervised_method == "IQR") {
      apply(methyl_data, 1L, IQR)
    } else {
      apply(methyl_data, 1L, mad)
    }

    rank_order <- order(metric, decreasing = TRUE)

    n_keep <- if (!is.null(top_n_features)) {
      message("Selecting top ", top_n_features, " features by ", unsupervised_method, ".")
      min(top_n_features, nrow(methyl_data))
    } else {
      message("Retaining top ", hvg_threshold, "% features by ", unsupervised_method, ".")
      ceiling((hvg_threshold / 100) * nrow(methyl_data))
    }

    filtered <- methyl_data[rank_order[seq_len(n_keep)], , drop = FALSE]

  } else if (unsupervised_method == "laplacian") {

    if (is.null(n_neighbors)) {
      n_neighbors <- max(2L, round(sqrt(ncol(methyl_data))))
      message("n_neighbors not set; using k = ", n_neighbors,
              " (approx. sqrt(n_samples)).")
    }

    message("Computing Laplacian scores (n_neighbors = ", n_neighbors, ").")

    ls_result  <- Rdimtools::do.lscore(
      X      = t(methyl_data),
      ndim   = nrow(methyl_data),
      nntype = "knn",
      params = list(k = n_neighbors)
    )

    ranked_idx <- ls_result$featidx   # ascending score = best first

    n_keep <- if (!is.null(top_n_features)) {
      message("Selecting top ", top_n_features, " features by Laplacian score.")
      min(top_n_features, length(ranked_idx))
    } else {
      message("Retaining top ", hvg_threshold, "% features by Laplacian score.")
      ceiling((hvg_threshold / 100) * nrow(methyl_data))
    }

    filtered <- methyl_data[ranked_idx[seq_len(n_keep)], , drop = FALSE]
  }

  message("Features after unsupervised filtering: ", nrow(filtered))
  filtered
}


# ══════════════════════════════════════════════════════════════════════════════
# CNV wrapper
# ══════════════════════════════════════════════════════════════════════════════

#' Prepare Copy Number Variation Data as SummarizedExperiment
#'
#' Processes CNV data by aligning samples with clinical data, removing constant
#' features, and optionally performing supervised, unsupervised, or biological
#' variable screening.
#'
#' @param cnv_data Numeric matrix of CNV data (genes/regions as rows, samples as
#'   columns).
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns). Must include \code{"overall_survival"} and
#'   \code{"status"} rows for supervised selection.
#' @param method Character string for variable screening: \code{"supervised"},
#'   \code{"unsupervised"}, \code{"biological"}, or \code{NULL} (no selection).
#'   Default is \code{NULL}.
#' @param supervised_method Character string specifying the supervised screening
#'   criterion, passed to \code{\link{supervised_screening}}:
#'   \itemize{
#'     \item \code{"pval"}: p-value of the marginal Cox test. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"padj"}: Benjamini-Hochberg adjusted p-value. Supports
#'       \code{p_threshold} and \code{top_n_features}.
#'     \item \code{"absmg"}: absolute marginal Cox coefficient
#'       \eqn{|\hat\beta_j|}. \strong{Requires \code{top_n_features}}.
#'     \item \code{"mgpval"}: combined score \eqn{-\log_{10}(p) \times
#'       |\hat\beta_j|}. \strong{Requires \code{top_n_features}}.
#'   }
#'   Default is \code{"pval"}.
#' @param unsupervised_method Character string specifying the unsupervised
#'   screening method: \code{"sparsity"} (fraction of zero entries) or
#'   \code{"laplacian"} (Laplacian score). Default is \code{"sparsity"}.
#' @param p_threshold Numeric p-value threshold for supervised selection. Used
#'   only when \code{supervised_method} is \code{"pval"} or \code{"padj"}
#'   and \code{top_n_features} is \code{NULL}. Default is 0.05.
#' @param max_zero_fraction Numeric proportion (0-1) for sparsity filtering:
#'   genes with a zero fraction above this value are removed. Used only when
#'   \code{unsupervised_method = "sparsity"} and \code{top_n_features} is
#'   \code{NULL}. Default is 0.75.
#' @param feature_fraction Numeric proportion (0-1) of features retained by
#'   Laplacian screening when \code{top_n_features} is \code{NULL}.
#'   Default is 0.25.
#' @param top_n_features Integer. If provided, overrides threshold/proportion
#'   filtering. \strong{Mandatory} when \code{supervised_method} is
#'   \code{"absmg"} or \code{"mgpval"}. Default is \code{NULL}.
#' @param gene_list Character vector of feature names to retain.
#'   Used only when \code{method = "biological"}. Default is \code{NULL}.
#' @param n_neighbors Integer. Number of neighbors for Laplacian score
#'   computation (\code{unsupervised_method = "laplacian"}). Defaults to
#'   \eqn{\approx \sqrt{n}} if \code{NULL}.
#' @param on_constant Character. After screening, any sample left constant
#'   (all-equal across the selected genes) is always repaired by adding back, for
#'   each such sample, the single original gene where it is most strongly altered,
#'   so that every retained sample has non-zero variance and downstream tools
#'   (e.g. flexynesis) do not drop it. Samples that are constant on the FULL input
#'   matrix carry no information and cannot be repaired; \code{on_constant}
#'   decides their fate: \code{"drop"} (default) removes them from the object,
#'   \code{"error"} stops and lists them, \code{"keep"} leaves them in (they stay
#'   constant).
#'
#' @return A \code{SummarizedExperiment} object with:
#'   \itemize{
#'     \item \code{assays$counts}: Filtered CNV matrix.
#'     \item \code{metadata$preprocessing}: List with \code{omics},
#'       \code{method}, \code{sub_method}, and number of features retained.
#'   }
#'
#' @export
prepare_CNV <- function(
    cnv_data,
    clinical_data,
    method              = NULL,
    supervised_method   = c("pval", "absmg", "padj", "mgpval"),
    unsupervised_method = c("sparsity", "laplacian"),
    p_threshold         = 0.05,
    max_zero_fraction   = 0.75,
    feature_fraction    = 0.25,
    top_n_features      = NULL,
    gene_list           = NULL,
    n_neighbors         = NULL,
    on_constant         = c("drop", "error", "keep")
) {

  supervised_method   <- match.arg(supervised_method)
  unsupervised_method <- match.arg(unsupervised_method)
  on_constant         <- match.arg(on_constant)
  if (!is.null(method)) {
    method <- match.arg(method, c("supervised", "unsupervised", "biological"))
  }

  samples_to_keep <- intersect(colnames(cnv_data), colnames(clinical_data))

  if (length(samples_to_keep) == 0) {
    stop("No overlapping samples between cnv_data and clinical_data.")
  }

  cnv_data      <- cnv_data[, samples_to_keep, drop = FALSE]
  clinical_data <- clinical_data[, samples_to_keep, drop = FALSE]

  message("Samples overlapping with clinical data: ", length(samples_to_keep))

  message("Removing constant features.")
  is_constant <- apply(cnv_data, 1L, function(x) length(unique(x)) == 1L)
  cnv_data    <- cnv_data[!is_constant, , drop = FALSE]

  message("Features remaining: ", nrow(cnv_data))

  # Keep the full (post-constant-feature-removal) matrix so the optional coverage
  # step can add back genes that screening discarded.
  cnv_full <- cnv_data

  if (is.null(method)) {
    message("No additional filtering applied.")

    return(
      SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = cnv_data),
        metadata = list(preprocessing = list(
          omics      = "CNV",
          method     = NULL,
          n_features = nrow(cnv_data)
        ))
      )
    )
  }

  cnv_data <- switch(
    method,
    supervised = {
      result <- supervised_screening(
        data           = cnv_data,
        clinical_data  = clinical_data,
        rank           = supervised_method,
        p_threshold    = p_threshold,
        top_n_features = top_n_features
      )
      result$data
    },

    unsupervised = .cnv_screen_unsupervised(
      cnv_data,
      unsupervised_method,
      max_zero_fraction,
      feature_fraction,
      top_n_features,
      n_neighbors
    ),

    biological = biological_screening(cnv_data, gene_list)$data,

    stop("Invalid method.")
  )

  # Always repair samples made constant by screening, so no patient is silently
  # dropped downstream.
  cov           <- .cnv_ensure_coverage(cnv_data, cnv_full, on_constant)
  cnv_data      <- cov$cnv
  clinical_data <- clinical_data[, colnames(cnv_data), drop = FALSE]

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = cnv_data),
    metadata = list(preprocessing = list(
      omics      = "CNV",
      method     = method,
      sub_method = if (method == "unsupervised")
        unsupervised_method else method,
      n_features = nrow(cnv_data)
    ))
  )

  return(se)
}


#' Number of features to keep for CNV screening
#'
#' Internal helper for \code{\link{prepare_CNV}}. Returns \code{top_n_features}
#' (capped at \code{n_total}) when supplied, otherwise a proportion of the total.
#'
#' @param n_total Integer. Total number of candidate features.
#' @param top_n_features Integer or \code{NULL}.
#' @param prop Numeric proportion (0-1) used when \code{top_n_features} is
#'   \code{NULL}.
#'
#' @return Integer number of features to retain.
#' @keywords internal
#' @noRd
.cnv_n_keep <- function(n_total, top_n_features, prop) {
  if (!is.null(top_n_features)) {
    min(as.integer(top_n_features), n_total)
  } else {
    max(1L, ceiling(prop * n_total))
  }
}


#' Unsupervised screening for CNV data
#'
#' Internal helper for \code{\link{prepare_CNV}}. Selects features by sparsity
#' (fraction of zero entries) or Laplacian score.
#'
#' @param cnv_data Numeric matrix (genes x samples).
#' @param unsupervised_method One of \code{"sparsity"}, \code{"laplacian"}.
#' @param max_zero_fraction Numeric proportion (0-1) for sparsity filtering.
#' @param feature_fraction Numeric proportion (0-1) retained by Laplacian
#'   screening when \code{top_n_features} is \code{NULL}.
#' @param top_n_features Integer or \code{NULL}. Overrides proportion-based
#'   filtering.
#' @param n_neighbors Integer or \code{NULL}. Neighbours for Laplacian score.
#'
#' @return Filtered numeric matrix (genes x samples).
#' @keywords internal
#' @noRd
.cnv_screen_unsupervised <- function(
    cnv_data,
    unsupervised_method,
    max_zero_fraction,
    feature_fraction,
    top_n_features,
    n_neighbors
) {

  message("Unsupervised screening method: ", unsupervised_method)
  n_total <- nrow(cnv_data)

  if (unsupervised_method == "sparsity") {

    zero_frac <- rowMeans(cnv_data == 0)

    if (!is.null(top_n_features)) {

      n_keep <- min(top_n_features, n_total)

      idx <- order(zero_frac, decreasing = FALSE)[seq_len(n_keep)]

      filtered <- cnv_data[idx, , drop = FALSE]

    } else {

      keep <- zero_frac <= max_zero_fraction

      message(
        "Kept ", sum(keep), " of ", n_total,
        " genes after sparsity filtering."
      )

      filtered <- cnv_data[keep, , drop = FALSE]
    }

  } else if (unsupervised_method == "laplacian") {

    if (is.null(n_neighbors)) {
      n_neighbors <- max(2L, round(sqrt(ncol(cnv_data))))
    }

    message("Computing Laplacian scores.")

    ls_result <- Rdimtools::do.lscore(
      X      = t(cnv_data),
      ndim   = nrow(cnv_data),
      nntype = "knn",
      params = list(k = n_neighbors)
    )

    ranked_idx <- ls_result$featidx

    n_keep <- .cnv_n_keep(
      length(ranked_idx),
      top_n_features,
      feature_fraction
    )

    filtered <- cnv_data[
      ranked_idx[seq_len(n_keep)],
      ,
      drop = FALSE
    ]
  }

  message("Features after unsupervised filtering: ", nrow(filtered))

  filtered
}


#' Repair samples left constant by CNV screening
#'
#' Internal helper for \code{\link{prepare_CNV}}. For each sample that becomes
#' constant across the \code{selected} genes, adds back the single gene from
#' \code{full} where that sample is most strongly altered. Samples constant
#' across \code{full} itself cannot be repaired and are handled per
#' \code{on_constant}.
#'
#' @param selected Screened matrix (selected genes x samples).
#' @param full Pre-screening matrix (all genes x same samples).
#' @param on_constant One of \code{"drop"}, \code{"error"}, \code{"keep"} for
#'   samples that are constant across the full matrix.
#'
#' @return A named list with \code{cnv} (repaired matrix), \code{added}
#'   (genes re-added), and \code{dropped} (uncoverable samples).
#' @keywords internal
#' @noRd
.cnv_ensure_coverage <- function(selected, full, on_constant) {

  is_flat <- function(m) apply(m, 2L, function(x) length(unique(x)) == 1L)

  # Samples uncoverable on the full matrix.
  full_flat <- colnames(full)[is_flat(full)]
  if (length(full_flat) > 0L) {
    if (on_constant == "error") {
      stop(length(full_flat), " sample(s) are constant across ALL genes and ",
           "cannot be made informative: ",
           paste(utils::head(full_flat, 15), collapse = ", "),
           if (length(full_flat) > 15)
             paste0(" ... (+", length(full_flat) - 15, " more)"),
           "\nSet on_constant = 'drop' to exclude them, or 'keep' to retain them.")
    } else if (on_constant == "drop") {
      keep_s   <- setdiff(colnames(selected), full_flat)
      selected <- selected[, keep_s, drop = FALSE]
      full     <- full[, keep_s, drop = FALSE]
      message("Coverage: dropped ", length(full_flat),
              " sample(s) constant across all genes.")
    }
    # "keep": leave them; they remain constant.
  }

  # Coverable samples left constant by screening: add their most-altered gene.
  flat_now <- colnames(selected)[is_flat(selected)]
  flat_now <- setdiff(flat_now, full_flat)
  flat_now <- intersect(flat_now, colnames(full))

  added <- character(0)
  for (s in flat_now) {
    g <- rownames(full)[which.max(abs(full[, s]))]
    if (!g %in% rownames(selected)) {
      selected <- rbind(selected, full[g, , drop = FALSE])
      rownames(selected)[nrow(selected)] <- g
      added <- union(added, g)
    }
  }

  if (length(flat_now) > 0L) {
    message("Coverage: repaired ", length(flat_now),
            " sample(s) constant after screening by adding ",
            length(added), " gene(s).")
  }

  message("Features after coverage re-addition: ", nrow(selected),
          " (", length(added), " added for sample coverage).")

  list(cnv = selected, added = added, dropped = full_flat)
}


# ══════════════════════════════════════════════════════════════════════════════
# Test-set alignment
# ══════════════════════════════════════════════════════════════════════════════

#' Prepare Test Omics Data as SummarizedExperiment
#'
#' Aligns a test omics matrix to the feature set of the corresponding training
#' \code{SummarizedExperiment}, retaining only the features selected during
#' training and ordering rows consistently. The result is
#' packaged as a \code{SummarizedExperiment} inheriting the preprocessing
#' metadata of the training object.
#'
#' @param test_data Numeric matrix of omics data for the test set (features as
#'   rows, samples as columns). Row names must be set.
#' @param train_se A \code{SummarizedExperiment} object produced by one of the
#'   \code{prepare_*} functions (mRNA, miRNA, CNV, Methylation). The features
#'   retained in \code{assays(train_se)$counts} define the target feature set.
#' @param on_constant Character. What to do with test samples that are constant
#'   across the (remaining) training features: \code{"drop"} (default) removes
#'   them, \code{"error"} stops, \code{"keep"} retains them.
#'
#' @return A \code{SummarizedExperiment} object with:
#'   \itemize{
#'     \item \code{assays$counts}: Test matrix restricted to training features,
#'      in the same row order as \code{train_se}, with
#'       constant samples handled per \code{on_constant}.
#'     \item \code{metadata$preprocessing}: Inherited from \code{train_se},
#'       with \code{split = "test"} and \code{n_features} updated.
#'   }
#'
#' @seealso \code{\link{prepare_RNA}}, \code{\link{prepare_methyl}},
#'   \code{\link{prepare_CNV}}
#'
#' @export
prepare_test_SE <- function(test_data, train_se,
                            on_constant = c("drop", "error", "keep")) {

  on_constant <- match.arg(on_constant)

  # ── Validate inputs ──────────────────────────────────────────────
  if (!inherits(train_se, "SummarizedExperiment")) {
    stop("'train_se' must be a SummarizedExperiment object.")
  }
  if (is.null(rownames(test_data))) {
    stop("'test_data' must have row names (feature identifiers).")
  }

  train_features <- rownames(SummarizedExperiment::assays(train_se)$counts)

  # ── Check for missing features ───────────────────────────────────
  missing_features <- setdiff(train_features, rownames(test_data))
  if (length(missing_features) > 0L) {
    stop(
      length(missing_features), " training feature(s) absent from test_data: ",
      paste(head(missing_features, 10L), collapse = ", "),
      if (length(missing_features) > 10L) " ..." else ""
    )
  }

  dropped_extra <- setdiff(rownames(test_data), train_features)
  if (length(dropped_extra) > 0L) {
    message(length(dropped_extra), " feature(s) in test_data not in training ",
            "set and will be dropped.")
  }

  # ── Align to training feature set ────────────────────────────────
  test_data <- test_data[train_features, , drop = FALSE]
  message("Test set aligned: ", nrow(test_data), " features x ",
          ncol(test_data), " samples.")

  # ── Handle samples constant across the training features ─────────
  layer <- tryCatch(
    S4Vectors::metadata(train_se)$preprocessing$omics,
    error = function(e) NULL
  )
  layer <- if (is.null(layer)) "this layer" else paste0("'", layer, "'")

  sd_per_sample <- apply(test_data, 2L, function(x) stats::sd(x, na.rm = TRUE))
  const         <- !is.finite(sd_per_sample) | sd_per_sample == 0

  if (any(const)) {
    bad <- colnames(test_data)[const]
    msg <- paste0(
      sum(const), " test sample(s) are constant across the training features ",
      "of ", layer, " (zero/NaN variance): ",
      paste(utils::head(bad, 15L), collapse = ", "),
      if (length(bad) > 15L) paste0(" ... (+", length(bad) - 15L, " more)") else ""
    )
    if (on_constant == "error") {
      stop(msg, "\nSet on_constant = 'drop' to exclude them, or 'keep' to retain them.")
    } else if (on_constant == "drop") {
      message(msg)
      test_data <- test_data[, !const, drop = FALSE]
      message("Dropped ", sum(const), " constant test sample(s); ",
              ncol(test_data), " remain.")
    } else {
      message(msg, " (kept; a downstream tool may drop them).")
    }
  }

  # ── Inherit and update metadata ──────────────────────────────────
  train_meta           <- S4Vectors::metadata(train_se)$preprocessing
  test_meta            <- train_meta
  test_meta$split      <- "test"
  test_meta$n_features <- nrow(test_data)

  # ── Build SummarizedExperiment ───────────────────────────────────
  se <- SummarizedExperiment::SummarizedExperiment(
    assays   = list(counts = test_data),
    metadata = list(preprocessing = test_meta)
  )

  return(se)
}
