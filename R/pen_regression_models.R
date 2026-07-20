#' Train Penalized Regression Survival Models on a MultiAssayExperiment Object
#'
#' Fits penalized Cox or AFT regression models on one or more omics layers from
#' a training \code{MultiAssayExperiment}, supporting single-omic, early
#' integration, and late integration strategies. Clinical data (survival outcome
#' and optional covariates) are extracted directly from \code{colData(mae_train)}.
#' Results are stored in the returned MAE.
#'
#' @param mae_train A \code{MultiAssayExperiment} object containing the training
#'   data. Must have \code{colData} with \code{"overall_survival"} and \code{"status"} columns.
#' @param omic_list Character vector of assay names in \code{mae_train} to use
#'   (e.g. \code{"mRNA"} or \code{c("mRNA", "CNV")}).
#' @param integration Character string specifying the integration strategy:
#'   \itemize{
#'     \item \code{"none"}: fits one model per omic in \code{omic_list}.
#'     \item \code{"early"}: concatenates all omics features into a single
#'       matrix and fits one model.
#'     \item \code{"late"}: fits one model per omic independently, then
#'       combines predictions via \code{late_method}.
#'   }
#'   Default is \code{"none"}.
#' @param model_method Character string specifying the penalized model:
#'   \code{"cox_lasso"}, \code{"cox_adaptive_lasso"},
#'   \code{"cox_elastic_net"}, or \code{"aft_pen"}. Default is
#'   \code{"cox_lasso"}.
#' @param late_method Character string for late integration ensemble:
#'   \code{"RRA"} (Robust Rank Aggregation) or \code{"average"}. Used only
#'   when \code{integration = "late"}. Default is \code{"RRA"}.
#' @param clinical_vars Character vector of column names in \code{colData(mae_train)}
#'   to append as additional covariates to the design matrix. Default is
#'   \code{NULL}.
#' @param seed Integer seed for reproducibility. Default is 123.
#' @param ... Additional arguments passed to the underlying train function
#'   (e.g. \code{nfolds}, \code{optimal_lambda}, \code{alpha}).
#'
#' @return \code{mae_train} with results stored in
#'   \code{metadata(mae_train)$analysis}:
#'   \itemize{
#'     \item \code{"none"} / \code{"early"}: per-assay results at
#'       \code{[[model_method]][[assay_name]]}, containing \code{train}
#'       (model, predictions, selected features) and \code{metrics}
#'       (resubstitution Harrell C-index).
#'     \item \code{"late"}: per-omic results as above, plus ensemble results
#'       at \code{$late_integration[[model_method]][[tag]]} (where \code{tag}
#'       encodes the integrated omics), containing \code{train},
#'       \code{metrics}, and \code{late_method}.
#'   }
#'
#' @seealso \code{\link{predict_pen_reg_models}}, \code{\link{build_MAE}},
#'   \code{\link{early_integration}}, \code{\link{late_integration}}
#'
#' @export
train_pen_reg_models <- function(
    mae_train,
    omic_list,
    integration   = c("none", "early", "late"),
    model_method  = c("cox_lasso", "cox_adaptive_lasso",
                      "cox_elastic_net", "aft_pen"),
    late_method   = c("RRA", "average"),
    clinical_vars = NULL,
    seed          = 123,
    ...
) {
  integration  <- match.arg(integration)
  model_method <- match.arg(model_method)
  late_method  <- match.arg(late_method)

  # ── Validate assay names ─────────────────────────────────────────
  missing_assays <- setdiff(omic_list, names(mae_train))
  if (length(missing_assays) > 0L) {
    stop("Assay(s) not found in mae_train: ",
         paste(missing_assays, collapse = ", "))
  }

  # ── Dispatch ─────────────────────────────────────────────────────
  if (integration == "early") {

    # nome univoco basato sulle omiche integrate
    tag <- paste0("early_", paste(sort(omic_list), collapse = "_"))

    mae_train <- early_integration(mae_train,
                                   omics_names    = omic_list,
                                   new_layer_name = tag)
    mae_train <- .fit_single(mae_train, tag,
                             clinical_vars, model_method, seed, ...)

  } else if (integration == "none") {

    for (assay_name in omic_list) {
      mae_train <- .fit_single(mae_train, assay_name,
                               clinical_vars, model_method, seed, ...)
    }

  } else if (integration == "late") {

    # ── fit per-omic ─────────────────────────────────────────────
    for (assay_name in omic_list) {
      already_fit <- !is.null(
        MultiAssayExperiment::metadata(mae_train)$analysis[[model_method]][[assay_name]]$train
      )
      if (already_fit) {
        message("Reusing existing fit for assay: ", assay_name)
      } else {
        message("Training assay: ", assay_name)
        mae_train <- .fit_single(
          mae_train, assay_name,
          clinical_vars, model_method, seed, ...
        )
      }
    }

    # ── collect per-omic train predictions ───────────────────────
    res_train_list <- lapply(omic_list, function(nm) {
      MultiAssayExperiment::metadata(mae_train)$analysis[[model_method]][[nm]]$train
    })
    names(res_train_list) <- omic_list

    # ── y on common patients (all same by construction) ───────────
    y_common <- .extract_xy(mae_train, omic_list[[1]], clinical_vars)$y

    # ── late integration ──────────────────────────────────────────
    late_train <- late_integration(
      res_list = res_train_list,
      y        = y_common,
      model    = if (model_method == "aft_pen") "aft_pen" else NULL,
      method   = late_method
    )

    surv_common <- survival::Surv(y_common$overall_survival, y_common$status)
    pred_late   <- late_train$pred

    metrics_train <- list(
      harrell_train = survival::concordance(
        surv_common ~ pred_late,
        reverse = TRUE
      )$concordance
    )

    # nome univoco basato sulle omiche integrate
    tag <- paste0("late_", paste(sort(omic_list), collapse = "_"))

    if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis$late_integration))
      MultiAssayExperiment::metadata(mae_train)$analysis$late_integration <- list()
    if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis$late_integration[[model_method]]))
      MultiAssayExperiment::metadata(mae_train)$analysis$late_integration[[model_method]] <- list()

    MultiAssayExperiment::metadata(mae_train)$analysis$late_integration[[model_method]][[tag]] <- list(
      train       = late_train,
      metrics     = metrics_train,
      late_method = late_method
    )
  }

  mae_train
}


#' Predict from Trained Penalized Regression Models on a New MAE
#'
#' Applies models trained by \code{\link{train_pen_reg_models}} to a new
#' \code{MultiAssayExperiment} (test set, validation cohort, or future data).
#' Survival outcome and clinical covariates are extracted directly from
#' \code{colData(mae_pred)} and \code{colData(mae_train)}.
#'
#' @param mae_pred A \code{MultiAssayExperiment} object containing the data to
#'   predict on, built with \code{\link{build_MAE}} after aligning features via
#'   \code{\link{prepare_test_SE}}.
#' @param mae_train A \code{MultiAssayExperiment} object containing fitted
#'   models in \code{metadata(mae_train)$analysis}, produced by
#'   \code{\link{train_pen_reg_models}}.
#' @param omic_list Character vector of assay names to predict on. Must match
#'   those used in \code{\link{train_pen_reg_models}}.
#' @param integration Character string specifying the integration strategy used
#'   during training: \code{"none"}, \code{"early"}, or \code{"late"}.
#'   Default is \code{"none"}.
#' @param model_method Character string specifying the model to retrieve:
#'   \code{"cox_lasso"}, \code{"cox_adaptive_lasso"},
#'   \code{"cox_elastic_net"}, or \code{"aft_pen"}. Must match the
#'   \code{model_method} used in training. Default is \code{"cox_lasso"}.
#' @param clinical_vars Character vector of column names in \code{colData}
#'   to append as covariates. Must match those used during training. Default
#'   is \code{NULL}.
#'
#' @return \code{mae_pred} with results stored in
#'   \code{metadata(mae_pred)$analysis}:
#'   \itemize{
#'     \item \code{"none"} / \code{"early"}: per-assay results at
#'       \code{[[model_method]][[assay_name]]}, containing \code{test}
#'       (predictions) and \code{metrics} (Harrell + Uno C-index).
#'     \item \code{"late"}: per-omic results as above, plus ensemble results
#'       at \code{$late_integration[[model_method]][[tag]]} (where \code{tag}
#'       encodes the integrated omics), containing \code{test} and
#'       \code{metrics}. The \code{late_method} is retrieved automatically
#'       from \code{mae_train} metadata.
#'   }
#'
#' @seealso \code{\link{train_pen_reg_models}}, \code{\link{prepare_test_SE}},
#'   \code{\link{build_MAE}}
#'
#' @export
predict_pen_reg_models <- function(
    mae_pred,
    mae_train,
    omic_list,
    integration  = c("none", "early", "late"),
    model_method = c("cox_lasso", "cox_adaptive_lasso",
                     "cox_elastic_net", "aft_pen"),
    clinical_vars = NULL
) {

  integration  <- match.arg(integration)
  model_method <- match.arg(model_method)

  # ── Validate models exist ────────────────────────────────────────
  if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis)) {
    stop("No fitted models found in mae_train metadata. ",
         "Run train_pen_reg_models() first.")
  }

  # ── Dispatch ─────────────────────────────────────────────────────
  if (integration == "early") {

    # stesso tag usato in training
    tag <- paste0("early_", paste(sort(omic_list), collapse = "_"))

    mae_pred <- early_integration(mae_pred,
                                  omics_names    = omic_list,
                                  new_layer_name = tag)
    mae_pred <- .predict_single(mae_pred, mae_train, tag,
                                clinical_vars, model_method)

  } else if (integration == "none") {

    for (assay_name in omic_list) {
      mae_pred <- .predict_single(mae_pred, mae_train, assay_name,
                                  clinical_vars, model_method)
    }

  } else if (integration == "late") {

    # ── predict per-omic ─────────────────────────────────────────
    for (assay_name in omic_list) {
      mae_pred <- .predict_single(mae_pred, mae_train, assay_name,
                                  clinical_vars, model_method)
    }

    # ── collect per-omic test predictions ────────────────────────
    res_test_list <- lapply(omic_list, function(nm) {
      MultiAssayExperiment::metadata(mae_pred)$analysis[[model_method]][[nm]]$test
    })
    names(res_test_list) <- omic_list

    # tag
    tag <- paste0("late_", paste(sort(omic_list), collapse = "_"))

    # ── late_method from mae_train (percorso con tag) ─────────────
    late_method <- MultiAssayExperiment::metadata(mae_train)$analysis$late_integration[[model_method]][[tag]]$late_method
    if (is.null(late_method)) {
      stop("Could not retrieve late_method from mae_train for '", tag, "'. ",
           "Ensure train_pen_reg_models() was run with integration = 'late' ",
           "on the same omic_list.")
    }

    # ── y_test and y_train from colData ───────────────────────────
    y_test  <- .extract_xy(mae_pred,  omic_list[[1]], clinical_vars = NULL)$y
    y_train <- .extract_xy(mae_train, omic_list[[1]], clinical_vars = NULL)$y

    # ── late integration on test ──────────────────────────────────
    late_test <- late_integration(
      res_list = res_test_list,
      y        = y_test,
      model    = if (model_method == "aft_pen") "aft_pen" else NULL,
      method   = late_method
    )

    surv_train <- survival::Surv(y_train$overall_survival, y_train$status)
    surv_test  <- survival::Surv(y_test$overall_survival,  y_test$status)
    pred_test  <- late_test$pred

    unotimes <- min(
      quantile(y_test$overall_survival, 0.90),
      max(y_test[y_test$status != 0, "overall_survival"])
    )

    metrics <- list(
      harrell_test = survival::concordance(
        surv_test ~ pred_test,
        reverse = TRUE
      )$concordance,
      uno_test = survAUC::UnoC(
        Surv.rsp     = surv_train,
        Surv.rsp.new = surv_test,
        lpnew        = pred_test,
        time         = unotimes
      )
    )

    # ── tag ──────────────────────────────
    if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis$late_integration))
      MultiAssayExperiment::metadata(mae_pred)$analysis$late_integration <- list()
    if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis$late_integration[[model_method]]))
      MultiAssayExperiment::metadata(mae_pred)$analysis$late_integration[[model_method]] <- list()

    MultiAssayExperiment::metadata(mae_pred)$analysis$late_integration[[model_method]][[tag]] <- list(
      test    = late_test,
      metrics = metrics
    )
  }

  mae_pred
}


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
  mae_matched <- MultiAssayExperiment::intersectColumns(mae[, , omics_names])

  # 2.
  list_of_matrices <- lapply(omics_names, function(name) {
    mat <- MultiAssayExperiment::assay(mae_matched[[name]])
    rownames(mat) <- paste0(name, "_", rownames(mat))
    return(mat)
  })

  combined_matrix <- do.call(rbind, list_of_matrices)

  # 3. SummarizedExperiment
  combined_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = combined_matrix))

  # 4.
  final_mae <- c(mae_matched,
                 setNames(list(combined_se), new_layer_name),
                 mapFrom = omics_names[1])

  return(final_mae)
}


#' Late Integration of Multi-Omics Survival Predictions
#'
#' Combines risk predictions from multiple independently-trained omics models into
#' a single consensus prediction. Supports averaging of normalized ranks or Robust
#' Rank Aggregation (RRA) to create an ensemble prediction that leverages information
#' from all omics layers.
#'
#' @param res_list A named list where each element contains results from a single
#'   omics layer model. Each element must include:
#'   \itemize{
#'     \item \code{pred}: Numeric vector of risk predictions for samples
#'   }
#'   List names typically represent omics layers (e.g., "transcriptomics", "proteomics").
#' @param y A data frame containing survival outcomes with columns:
#'   \itemize{
#'     \item \code{overall_survival}: Numeric survival time
#'     \item \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#'   Row names should correspond to patient/sample identifiers.
#' @param model Character string indicating whether the underlying prediction
#'   model is a Penalized Accelerated Failure Time model (\code{"aft_pen"}).
#'   Default is \code{NULL}, corresponding to the use of Cox models.
#' @param method Character string specifying the integration method. Must be one of:
#'   \itemize{
#'     \item \code{"average"}: Average of normalized ranks across omics layers
#'     \item \code{"RRA"}: Robust Rank Aggregation using the RobustRankAggreg package
#'   }
#'
#' @return A list containing:
#'   \item{pred}{Numeric vector of consensus risk ranks for all patients. Lower ranks
#'     indicate lower risk (longer predicted survival). Length equals number of patients.}
#'   \item{merged_risks}{Data frame containing:
#'     \itemize{
#'       \item \code{patient_id}: Sample identifiers
#'       \item \code{risk_omics1}, \code{risk_omics2}, ...: Original risk scores from each omics layer
#'       \item \code{rank_omics1}, \code{rank_omics2}, ...: Ranks of normalized scores within each layer
#'       \item \code{consensus_rank}: Final integrated consensus rank
#'     }}
#'
#' @details
#' The late integration workflow consists of the following steps:
#' \enumerate{
#'   \item \strong{Extract}: Collect risk scores from each omics layer model
#'   \item \strong{Merge}: Combine risk scores into a single data frame by patient ID
#'   \item \strong{Rank}: Convert normalized scores to ranks within each layer using
#'     average method for ties
#'   \item \strong{Integrate}: Combine ranks across layers:
#'     \itemize{
#'       \item \code{average}: Simple arithmetic mean of ranks
#'       \item \code{RRA}: Statistical aggregation accounting for rank distribution
#'     }
#' }
#'
#' @seealso
#' \code{\link[RobustRankAggreg]{aggregateRanks}} for RRA implementation
#'
#' @export
late_integration <- function(res_list, y,
                             model  = NULL,
                             method = c("average", "RRA")) {

  method <- match.arg(method)

  is_aft      <- !is.null(model) && model == "aft_pen"
  patient_ids <- rownames(y)
  n           <- nrow(y)

  # 1. Collect risk scores per omic
  risk_dfs <- lapply(seq_along(res_list), function(i) {
    df <- data.frame(patient_id = patient_ids,
                     risk_score = as.numeric(res_list[[i]]$pred))
    colnames(df)[2] <- paste0("risk_omics", i)
    df
  })

  # 2. Merge + re-align to original patient order
  merged_risks <- Reduce(function(a, b) merge(a, b, by = "patient_id"),
                         risk_dfs)
  merged_risks <- merged_risks[match(patient_ids, merged_risks$patient_id), ]

  # 3. Rank each omic (direction-aware)
  rank_cols <- paste0("rank_omics", seq_along(res_list))

  for (i in seq_along(res_list)) {
    raw <- merged_risks[[paste0("risk_omics", i)]]
    # AFT: higher LP = longer survival = lower risk -> negate before ranking
    # Cox: higher LP = higher risk -> rank as-is
    merged_risks[[rank_cols[i]]] <- rank(
      if (is_aft) -raw else raw,
      ties.method = "average"
    )
  }

  # 4. Consensus
  if (method == "average") {

    merged_risks$consensus_rank <- rowMeans(
      merged_risks[, rank_cols, drop = FALSE]
    )

  } else if (method == "RRA") {

    rank_matrix <- as.matrix(merged_risks[, rank_cols, drop = FALSE])
    rownames(rank_matrix) <- merged_risks$patient_id

    rra_result <- RobustRankAggreg::aggregateRanks(
      rmat = rank_matrix / n,
      N    = n
    )

    # re-align RRA scores to original patient order
    ord <- match(merged_risks$patient_id, rra_result$Name)
    merged_risks$consensus_score <- rra_result$Score[ord]

    # rank 1 = lowest aggregated score = lowest risk
    merged_risks$consensus_rank <- rank(
      merged_risks$consensus_score, ties.method = "average"
    )
  }

  return(list(
    pred         = merged_risks$consensus_rank,   # rank 1 = lowest risk
    merged_risks = merged_risks
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Fit a single penalized model on one assay and store it in the MAE
#'
#' Internal helper for \code{\link{train_pen_reg_models}}. Extracts the design
#' matrix and survival outcome for a single assay, dispatches to the appropriate
#' \code{train_*} function, computes the resubstitution Harrell C-index, and
#' stores the fit and metrics in \code{metadata(mae_train)$analysis}.
#'
#' @param mae_train A \code{MultiAssayExperiment} object.
#' @param assay_name Character. Name of the assay to fit.
#' @param clinical_vars Character vector of clinical covariates, or \code{NULL}.
#' @param model_method Character. One of \code{"cox_lasso"},
#'   \code{"cox_adaptive_lasso"}, \code{"cox_elastic_net"}, \code{"aft_pen"}.
#' @param seed Integer seed for reproducibility.
#' @param ... Additional arguments passed to the underlying \code{train_*} function.
#'
#' @return The input \code{mae_train} with the fit stored in metadata.
#' @keywords internal
#' @noRd
.fit_single <- function(mae_train, assay_name,
                        clinical_vars, model_method, seed, ...) {

  train_data <- .extract_xy(mae_train, assay_name, clinical_vars)

  train_fun <- switch(model_method,
                      cox_lasso          = train_cox_glmnet_LASSO,
                      cox_adaptive_lasso = train_cox_glmnet_adaptiveLASSO,
                      cox_elastic_net    = train_cox_glmnet_ElasticNet,
                      aft_pen            = train_penAFT
  )

  train_res  <- train_fun(train_data$x, train_data$y, seed = seed, ...)
  pred_train <- if (model_method == "aft_pen") -train_res$pred else train_res$pred
  surv_train <- survival::Surv(train_data$y$overall_survival,
                               train_data$y$status)

  metrics <- list(
    harrell_train = survival::concordance(
      surv_train ~ pred_train,
      reverse = TRUE
    )$concordance
  )

  if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis))
    MultiAssayExperiment::metadata(mae_train)$analysis <- list()
  if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis[[model_method]]))
    MultiAssayExperiment::metadata(mae_train)$analysis[[model_method]] <- list()

  MultiAssayExperiment::metadata(mae_train)$analysis[[model_method]][[assay_name]] <- list(
    train   = train_res,
    metrics = metrics
  )

  mae_train
}


#' Extract design matrix and survival outcome for one assay
#'
#' Internal helper shared by \code{\link{train_pen_reg_models}} and
#' \code{\link{predict_pen_reg_models}}. Transposes the assay counts to a
#' samples-by-features matrix, pulls \code{"overall_survival"} and
#' \code{"status"} from \code{colData}, and optionally appends clinical
#' covariates (factors coerced to numeric) as leading columns of the design
#' matrix.
#'
#' @param mae A \code{MultiAssayExperiment} object.
#' @param assay_name Character. Name of the assay to extract.
#' @param clinical_vars Character vector of clinical covariates to append, or
#'   \code{NULL} (default).
#'
#' @return A named list with \code{x} (design matrix, samples in rows) and
#'   \code{y} (data frame with \code{overall_survival} and \code{status}).
#' @keywords internal
#' @noRd
.extract_xy <- function(mae, assay_name, clinical_vars = NULL) {

  x       <- t(SummarizedExperiment::assays(mae[[assay_name]])$counts)
  samples <- rownames(x)

  clin <- as.data.frame(
    MultiAssayExperiment::colData(mae)[samples, , drop = FALSE]
  )

  # ── Validate survival columns ────────────────────────────────────
  required_cols <- c("overall_survival", "status")
  missing_cols  <- setdiff(required_cols, colnames(clin))
  if (length(missing_cols) > 0L) {
    stop("colData(mae) missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  y <- data.frame(
    overall_survival = as.numeric(clin$overall_survival),
    status           = as.numeric(clin$status),
    row.names        = samples
  )

  if (!is.null(clinical_vars)) {
    missing_vars <- setdiff(clinical_vars, colnames(clin))
    if (length(missing_vars) > 0L) {
      stop("Clinical variable(s) not found in colData(mae): ",
           paste(missing_vars, collapse = ", "))
    }
    clin_subset <- clin[, clinical_vars, drop = FALSE]
    clin_subset <- as.data.frame(lapply(clin_subset, function(v) {
      if (is.factor(v)) as.numeric(v) else v
    }))
    x <- as.matrix(cbind(clin_subset, x))
  }

  list(x = x, y = y)
}


#' Predict on a single assay and store results in the MAE
#'
#' Internal helper for \code{\link{predict_pen_reg_models}}. Retrieves the
#' fitted model for one assay from \code{mae_train}, applies the appropriate
#' \code{test_*} function to \code{mae_pred}, computes Harrell and Uno
#' C-indices, and stores predictions and metrics in
#' \code{metadata(mae_pred)$analysis}.
#'
#' @param mae_pred A \code{MultiAssayExperiment} object to predict on.
#' @param mae_train A \code{MultiAssayExperiment} object holding the fitted models.
#' @param assay_name Character. Name of the assay to predict.
#' @param clinical_vars Character vector of clinical covariates, or \code{NULL}.
#' @param model_method Character. One of \code{"cox_lasso"},
#'   \code{"cox_adaptive_lasso"}, \code{"cox_elastic_net"}, \code{"aft_pen"}.
#'
#' @return The input \code{mae_pred} with predictions stored in metadata.
#' @keywords internal
#' @noRd
.predict_single <- function(mae_pred, mae_train, assay_name,
                            clinical_vars, model_method) {

  # ── retrieve fitted model from mae_train ──────────────────────────
  train_res <- MultiAssayExperiment::metadata(mae_train)$analysis[[model_method]][[assay_name]]$train
  if (is.null(train_res)) {
    stop("No fitted model found for assay '", assay_name,
         "' and method '", model_method, "' in mae_train.")
  }

  # ── y_train and y_test from colData ───────────────────────────────
  y_train <- .extract_xy(mae_train, assay_name, clinical_vars = NULL)$y
  test_data <- .extract_xy(mae_pred, assay_name, clinical_vars)

  test_fun <- switch(model_method,
                     cox_lasso          = ,
                     cox_adaptive_lasso = ,
                     cox_elastic_net    = test_cox_glmnet,
                     aft_pen            = test_penAFT
  )

  test_res   <- test_fun(test_data$x, train_res)
  pred_test <- if (model_method == "aft_pen") -test_res$pred else test_res$pred
  surv_train <- survival::Surv(y_train$overall_survival, y_train$status)
  surv_test  <- survival::Surv(test_data$y$overall_survival, test_data$y$status)

  unotimes <- min(
    quantile(test_data$y$overall_survival, 0.90),
    max(test_data$y[test_data$y$status != 0, "overall_survival"])
  )

  metrics <- list(
    harrell_test = survival::concordance(
      surv_test ~ pred_test,
      reverse = TRUE
    )$concordance,
    uno_test = survAUC::UnoC(
      Surv.rsp     = surv_train,
      Surv.rsp.new = surv_test,
      lpnew        = pred_test,
      time         = unotimes
    )
  )

  if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis))
    MultiAssayExperiment::metadata(mae_pred)$analysis <- list()
  if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis[[model_method]]))
    MultiAssayExperiment::metadata(mae_pred)$analysis[[model_method]] <- list()

  MultiAssayExperiment::metadata(mae_pred)$analysis[[model_method]][[assay_name]] <- list(
    test    = test_res,
    metrics = metrics
  )

  mae_pred
}
