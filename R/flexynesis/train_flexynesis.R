train_flexynesis <- function(
    X_train_list,
    X_test_list,
    y_train,
    y_test,
    omic_list,
    surv_time = "overall_survival",
    surv_event = "status",
    tmpdir = tempdir(),
    conda_path = NULL,
    ...
) {

  configure_flexynesis_python(conda_path = conda_path)

  export_flexynesis_data(
    X_train_list, X_test_list,
    y_train, y_test,
    outdir = tmpdir
  )

  run_flexynesis(
    data_path = tmpdir,
    omic_list = omic_list,
    surv_time = surv_time,
    surv_event = surv_event,
    outdir = file.path(tmpdir, "results"),
    ...
  )

  preds <- import_flexynesis_results(
    file.path(tmpdir, "results", "job.predicted_labels.csv")
  )

  train_pred <- preds$train
  test_pred <- preds$test
  surv_obj <- survival::Surv(y_train[[surv_time]], y_train[[surv_event]])

  cindex <- survival::concordance(
    surv_obj ~ train_pred$predicted_label,
    reverse = TRUE
  )$concordance

  list(
    model = "flexynesis",
    pred = train_pred,
    pred_test= test_pred,
    cindex = cindex,
    selected_features = NULL,
    params = list(...)
  )
}
