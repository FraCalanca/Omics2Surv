plot_risk_scores = function (test_list, color="green", bins= 40, alpha = 0.6, title="Risk scores distribution", data_type= NULL) {

  if (!"pred_test" %in% names(test_list)) {
    stop("Error: pred_test is not in test_list")
  }

  if (!is.null(data_type)){
    title = paste(title, " - ", data_type)
  }

  ggplot2::ggplot(data.frame(score = as.numeric(test_list$pred_test)), ggplot2::aes(x = score)) +
    ggplot2::geom_histogram(bins = bins, fill = color, alpha = alpha, color = "black") +
    ggplot2::labs(title = title, x = "Risk Score", y = "Count") +
    ggplot2::theme_minimal()
}
