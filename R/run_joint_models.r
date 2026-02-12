#' Run Joint Multi-Omics Survival Models
#'
#' Performs end-to-end joint modeling workflow: builds a MultiAssayExperiment from
#' omics data, prepares design matrices, splits data into train/test sets, and fits
#' a specified joint integration model.
#'
#' @param omics_list Named list of SummarizedExperiment objects, each representing
#'   an omics layer (e.g., transcriptomics, proteomics).
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns).
#' @param surv_list Character vector of survival variable names in clinical data
#'   (e.g., \code{c("overall_survival", "status")}).
#' @param joint_method Character string specifying the joint model: \code{"CoopCox"}
#'   (2-3 omics), \code{"AFTCoop"} (exactly 2 omics), \code{"blockForest"} (any number),
#'   or \code{"flexynesis"}.
#' @param train_fraction Numeric proportion (0-1) of samples for training. Default is 0.7.
#' @param seed Integer seed for reproducibility. Default is 123.
#' @param ... Additional arguments passed to the specific joint modeling function.
#'
#' @return A MultiAssayExperiment object with results stored in
#'   \code{metadata(mae)$analysis$cooperative[[joint_method]]} containing:
#'   \itemize{
#'     \item \code{train}: Training results (model, predictions, C-index, selected features)
#'     \item \code{test}: Test results (predictions, C-index)
#'   }
#'
#' @export
run_joint_models <- function(
    omics_list,
    clinical_data,
    surv_list,
    joint_method = c("CoopCox", "AFTCoop", "blockForest", "flexynesis"),
    train_fraction = 0.7,
    seed = 123,
    ...
) {

  joint_method <- match.arg(joint_method)

  # 1. Build MAE
  Mae <- build_MAE(omics_list, clinical_data)

  # 2. Design matrices
  X_list <- lapply(names(omics_list), function(omic) {
    prepare_design_matrix(Mae, assay_name = omic)
  })
  names(X_list) <- names(omics_list)

  # 3. Keep common patients
  common_ids <- Reduce(intersect, lapply(X_list, rownames))
  X_list <- lapply(X_list, function(x) x[common_ids, , drop = FALSE])

  # 4. Survival
  y <- prepare_surv_data(Mae, surv_list, sample_order = common_ids)

  # 5. Split
  split <- split_train_test(
    X_list[[1]],
    y,
    train_fraction = train_fraction,
    seed = seed
  )

  # 6. Add split info for joint models
  split$X_list <- X_list

  # 7. Fit joint model
  Mae <- J_models(
    mae = Mae,
    split = split,
    joint_method = joint_method,
    seed = seed,
    omics_list,
    ...
  )

  return(Mae)
}
