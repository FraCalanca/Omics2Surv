#' Train Joint Multi-Omics Survival Models
#'
#' Fits joint multi-omics survival models (\code{CoopCox}, \code{AFTCoop}, or
#' \code{blockForest}) on two or more omics layers from a training
#' \code{MultiAssayExperiment}. Survival outcome (\code{overall_survival} and
#' \code{status}) is extracted directly from \code{colData(mae_train)}. Results
#' are stored in the returned MAE.
#'
#' @param mae_train A \code{MultiAssayExperiment} built with
#'   \code{\link{build_MAE}}. \code{colData} must contain
#'   \code{"overall_survival"} and \code{"status"}.
#' @param omic_list Character vector of assay names to use
#'   (e.g. \code{c("mRNA", "CNV")}). Order matters for \code{"AFTCoop"}
#'   (exactly 2 omics required).
#' @param joint_method Character string: \code{"CoopCox"} (2-3 omics),
#'   \code{"AFTCoop"} (exactly 2 omics), or \code{"blockForest"} (any number).
#' @param seed Integer seed. Default is 123.
#' @param ... Additional arguments passed to the underlying train function:
#'   \itemize{
#'     \item \code{CoopCox}: \code{val_fraction}, \code{lambdaVector},
#'       \code{folds}, \code{alpha}, \code{maxit}.
#'     \item \code{AFTCoop}: \code{model}, \code{rho_values},
#'       \code{ncore_max_rho}, \code{ncore_max_cv}.
#'     \item \code{blockForest}: \code{threads},
#'       \code{always.split.variables}, \code{n_top_features}.
#'   }
#'
#' @return \code{mae_train} with results stored at
#'   \code{metadata(mae_train)$analysis$joint[[joint_method]]}, containing
#'   \code{train} and \code{metrics} (resubstitution Harrell C-index).
#'
#' @seealso \code{\link{predict_j_models}}, \code{\link{train_cooplearning}},
#'   \code{\link{train_AFTCoop}}, \code{\link{train_blockForest}}
#'
#' @export
train_j_models <- function(
    mae_train,
    omic_list,
    joint_method = c("CoopCox", "AFTCoop", "blockForest"),
    seed         = 123,
    ...
) {

  joint_method <- match.arg(joint_method)

  # ── Validate assay names ─────────────────────────────────────────
  missing_assays <- setdiff(omic_list, names(mae_train))
  if (length(missing_assays) > 0L) {
    stop("Assay(s) not found in mae_train: ",
         paste(missing_assays, collapse = ", "))
  }

  # ── Method-specific input checks ─────────────────────────────────
  if (joint_method == "AFTCoop" && length(omic_list) != 2L) {
    stop("'AFTCoop' requires exactly 2 omics in omic_list.")
  }
  if (joint_method == "CoopCox" && !length(omic_list) %in% c(2L, 3L)) {
    stop("'CoopCox' requires 2 or 3 omics in omic_list.")
  }

  # ── Extract from MAE ─────────────────────────────────────────────
  X_list <- .extract_x_list(mae_train, omic_list)
  y      <- .extract_xy(mae_train, omic_list[[1]])$y

  # ── Dispatch ─────────────────────────────────────────────────────
  train_res <- switch(
    joint_method,

    CoopCox = train_cooplearning(
      X_list     = X_list,
      y_training = y,
      seed       = seed,
      ...
    ),

    AFTCoop = train_AFTCoop(
      x_training1 = X_list[[1]],
      x_training2 = X_list[[2]],
      y_training  = y,
      seed        = seed,
      ...
    ),

    blockForest = train_blockForest(
      X_list     = X_list,
      y_training = y,
      seed       = seed,
      ...
    )
  )

  # ── Metrics ───────────────────────────────────────────────────────
  metrics <- .compute_train_metrics(train_res, y, joint_method)

  # ── Store in MAE metadata ─────────────────────────────────────────
  tag <- paste(sort(omic_list), collapse = "_")

  if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis))
    MultiAssayExperiment::metadata(mae_train)$analysis <- list()
  if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis$joint))
    MultiAssayExperiment::metadata(mae_train)$analysis$joint <- list()
  if (is.null(MultiAssayExperiment::metadata(mae_train)$analysis$joint[[joint_method]]))
    MultiAssayExperiment::metadata(mae_train)$analysis$joint[[joint_method]] <- list()

  MultiAssayExperiment::metadata(mae_train)$analysis$joint[[joint_method]][[tag]] <- list(
    train   = train_res,
    metrics = metrics
  )

  mae_train
}


#' Predict from Trained Joint Multi-Omics Models
#'
#' Applies joint models trained by \code{\link{train_j_models}} to a new MAE.
#' The underlying test functions (\code{\link{test_cooplearning}},
#' \code{\link{test_AFTCoop}}, \code{\link{test_blockForest}}) can always be
#' called directly.
#'
#' @param mae_pred A \code{MultiAssayExperiment} to predict on.
#' @param mae_train A \code{MultiAssayExperiment} with fitted models in
#'   \code{metadata(mae_train)$analysis$joint}.
#' @param omic_list Character vector of assay names. Must match those used in
#'   \code{\link{train_j_models}}.
#' @param joint_method Character string: \code{"CoopCox"}, \code{"AFTCoop"},
#'   or \code{"blockForest"}.
#' @param ... Additional arguments passed to the test function (e.g.
#'   \code{threads} for blockForest).
#'
#' @return \code{mae_pred} with results stored at
#'   \code{metadata(mae_pred)$analysis$joint[[joint_method]]}, containing
#'   \code{test} and \code{metrics} (Harrell + Uno C-index).
#'
#' @seealso \code{\link{train_j_models}}, \code{\link{test_cooplearning}},
#'   \code{\link{test_AFTCoop}}, \code{\link{test_blockForest}}
#'
#' @export
predict_j_models <- function(
    mae_pred,
    mae_train,
    omic_list,
    joint_method = c("CoopCox", "AFTCoop", "blockForest"),
    ...
) {

  joint_method <- match.arg(joint_method)

  # ── Retrieve fitted model ─────────────────────────────────────────
  tag <- paste(sort(omic_list), collapse = "_")

  train_res <- MultiAssayExperiment::metadata(mae_train)$analysis$joint[[joint_method]][[tag]]$train
  if (is.null(train_res)) {
    stop("No fitted '", joint_method, "' model found in mae_train for omics '",
         tag, "'. Run train_j_models() with the same omic_list first.")
  }

  # ── Extract from MAE ─────────────────────────────────────────────
  X_list  <- .extract_x_list(mae_pred, omic_list)
  y_test  <- .extract_xy(mae_pred,  omic_list[[1]])$y
  y_train <- .extract_xy(mae_train, omic_list[[1]])$y

  # ── Dispatch ─────────────────────────────────────────────────────
  test_res <- switch(
    joint_method,

    CoopCox = test_cooplearning(
      X_list_test    = X_list,
      training_result = train_res
    ),

    AFTCoop = test_AFTCoop(
      x_test1        = X_list[[1]],
      x_test2        = X_list[[2]],
      training_result = train_res
    ),

    blockForest = test_blockForest(
      X_list         = X_list,
      training_result = train_res,
      ...
    )
  )

  # ── Metrics ───────────────────────────────────────────────────────
  metrics <- .compute_test_metrics(train_res, test_res, y_train, y_test,
                                   joint_method)

  # ── Store in MAE metadata ─────────────────────────────────────────
  if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis))
    MultiAssayExperiment::metadata(mae_pred)$analysis <- list()
  if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis$joint))
    MultiAssayExperiment::metadata(mae_pred)$analysis$joint <- list()
  if (is.null(MultiAssayExperiment::metadata(mae_pred)$analysis$joint[[joint_method]]))
    MultiAssayExperiment::metadata(mae_pred)$analysis$joint[[joint_method]] <- list()

  MultiAssayExperiment::metadata(mae_pred)$analysis$joint[[joint_method]][[tag]] <- list(
    test    = test_res,
    metrics = metrics
  )

  mae_pred
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

.extract_x_list <- function(mae, omic_list) {
  lapply(setNames(omic_list, omic_list), function(nm) {
    t(SummarizedExperiment::assays(mae[[nm]])$counts)
  })
}

# Train: only harrell_train
# AFTCoop: pred is a list (one per rho) → compute per rho
.compute_train_metrics <- function(train_res, y, joint_method) {

  is_aft  <- joint_method == "AFTCoop"
  surv_y  <- survival::Surv(y$overall_survival, y$status)

  .harrell_train <- function(pred) {
    p <- if (is_aft) -pred else pred
    list(harrell_train = survival::concordance(surv_y ~ p, reverse = TRUE)$concordance)
  }

  if (is_aft) {
    lapply(train_res$pred, .harrell_train)
  } else {
    .harrell_train(train_res$pred)
  }
}

# Test: harrell_test + uno_test
# AFTCoop: one set of metrics per rho
.compute_test_metrics <- function(train_res, test_res, y_train, y_test,
                                  joint_method) {

  is_aft <- joint_method == "AFTCoop"

  .metrics_pair <- function(pred_train, pred_test) {
    .test_metrics(y_train, y_test, pred_train, pred_test, is_aft)
  }

  if (is_aft) {
    rho_names <- names(test_res$pred)
    setNames(
      lapply(rho_names, function(rho_name) {
        .metrics_pair(train_res$pred[[rho_name]], test_res$pred[[rho_name]])
      }),
      rho_names
    )
  } else {
    .metrics_pair(train_res$pred, test_res$pred)
  }
}

# Harrell test + Uno test
# Pass RAW predictions: AFT negation is
# handled internally via `is_aft`.
.test_metrics <- function(y_train, y_test, pred_train, pred_test, is_aft) {

  # direction: AFT higher LP = longer survival = lower risk -> negate
  if (is_aft) {
    pred_train <- -pred_train
    pred_test  <- -pred_test
  }

  surv_train <- survival::Surv(y_train$overall_survival, y_train$status)
  surv_test  <- survival::Surv(y_test$overall_survival,  y_test$status)

  harrell_test <- survival::concordance(
    surv_test ~ pred_test,
    reverse = TRUE
  )$concordance

  # Uno C: needs events in both sets; truncation time bounded by the train
  # follow-up (IPCW estimated on train) and the test event support.
  ev_test  <- y_test$overall_survival[y_test$status != 0]
  ev_train <- y_train$overall_survival[y_train$status != 0]

  if (length(ev_test) == 0L || length(ev_train) == 0L) {
    uno_test <- NA_real_
  } else {
    uno_time <- min(
      stats::quantile(y_test$overall_survival, 0.90, names = FALSE),
      max(ev_test),
      max(ev_train)
    )
    uno_test <- survAUC::UnoC(
      Surv.rsp     = surv_train,
      Surv.rsp.new = surv_test,
      lpnew        = as.numeric(pred_test),
      time         = uno_time
    )
  }

  list(
    harrell_test = harrell_test,
    uno_test     = uno_test
  )
}
