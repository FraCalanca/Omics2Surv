test_penAFT = function(x_test, y_test, training_result) {

  X_test <- as.matrix(x_test)
  time <- as.numeric(y_test$overall_survival)
  event <- as.numeric(y_test$status)

  # --- penAFT.predict ---
  linear_pred <- as.vector(
    penAFT::penAFT.predict(
      training_result$best_cv_fit,
      Xnew = X_test,
    )
  )

  # ---C-index---
  cindex_test= survival::concordance(survival::Surv(time, event) ~ linear_pred)

  # --- Output ---
  return(list(
    pred = linear_pred,
    cindex = cindex_test$concordance
  ))
}
