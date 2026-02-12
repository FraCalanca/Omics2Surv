#' Dichotomize Continuous Risk Scores into High and Low Risk Groups
#'
#' Converts continuous risk predictions into binary risk groups (High/Low) based on
#' a cutoff threshold. The cutoff can be determined from training data using median,
#' or a custom value, and is then applied to test data for risk stratification.
#'
#' @param pred_train Numeric vector of risk predictions from the training set.
#'   Used to determine the cutoff threshold. NA values are removed when calculating
#'   median cutoffs.
#' @param pred_test Numeric vector of risk predictions from the test set (or any
#'   dataset to be stratified). These predictions will be classified as High or Low
#'   risk based on the cutoff derived from \code{pred_train}.
#' @param method Character string specifying the method to determine the cutoff threshold.
#'   Options are:
#'   \itemize{
#'     \item \code{"median"}: Use the median of training predictions
#'     \item \code{"custom"}: Use a user-specified cutoff value (requires \code{custom_cutoff})
#'   }
#'   Default is \code{"median"}.
#' @param custom_cutoff Numeric value specifying a custom cutoff threshold. Only used
#'   when \code{method = "custom"}. Must be provided if custom method is selected.
#'   Default is NULL.
#'
#' @return A list containing:
#'   \item{groups}{Factor vector with two levels (\code{"Low"}, \code{"High"}) indicating
#'     risk group assignment for each sample in \code{pred_test}. Samples with predictions
#'     above the cutoff are classified as "High" risk; those at or below are "Low" risk.}
#'   \item{cutoff}{Numeric value of the threshold used for dichotomization. This allows
#'     users to see the exact cutoff applied and use it consistently across multiple datasets.}
#'
#' @details
#' This function is commonly used in survival analysis to:
#' \itemize{
#'   \item Stratify patients into risk groups for Kaplan-Meier analysis
#'   \item Create interpretable clinical risk categories from continuous scores
#'   \item Enable comparison of survival curves between high and low risk groups
#'   \item Facilitate risk-based treatment decision making
#' }
#'
#' Classification logic:
#' \itemize{
#'   \item \code{pred_test > cutoff} → "High" risk
#'   \item \code{pred_test <= cutoff} → "Low" risk
#' }
#'
#' @section Cutoff Methods:
#' \describe{
#'   \item{Median}{Robust to outliers; splits training samples into equal-sized groups.}
#'   \item{Custom}{Allows domain-specific thresholds based on clinical relevance,
#'     prior studies, or optimization criteria.}
#' }
#'
#' @export
dichotomize_risk_groups = function(pred_train, pred_test, method = "median", custom_cutoff = NULL) {

  cutoff <- switch(method,
                   median = median(pred_train, na.rm = TRUE),
                   mean = mean(pred_train, na.rm = TRUE),
                   custom = custom_cutoff,
                   stop("Invalid method. Choose 'median', 'mean', or 'custom'."))


  group <- ifelse(pred_test > cutoff, "High", "Low")

  return(list(groups=factor(group, levels = c("Low", "High")),cutoff=cutoff))
}
