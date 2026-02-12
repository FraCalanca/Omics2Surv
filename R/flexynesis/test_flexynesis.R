test_flexynesis <- function(
    training_result,
    y_test,
    surv_time = "overall_survival",
    surv_event = "status"
) {
  test_pred <- training_result$pred_test

  surv_obj <- survival::Surv(y_test[[surv_time]], y_test[[surv_event]])

  cindex <- survival::concordance(
    surv_obj ~ test_pred$predicted_label,
    reverse = TRUE
  )$concordance

  list(
    pred = test_pred,
    cindex = cindex
  )
}
