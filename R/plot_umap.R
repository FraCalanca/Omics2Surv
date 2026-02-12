plot_umap = function(data,
                                features,
                                pred_train = NULL,
                                pred_test = NULL,
                                cutoff = NULL,
                                color_by = c("none", "risk"),
                                palette = c("High" = "tomato", "Low" = "deepskyblue"),
                                base_color = "lightgreen") {
  # Controllo pacchetti
  if (!requireNamespace("umap", quietly = TRUE)) stop("Install 'umap' package.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2' package.")

  color_by = match.arg(color_by)

  # Scaling + feature selection
  data_scaled = scale(data)
  data_selected = data_scaled[, features, drop = FALSE]

  # UMAP
  umap_res = umap(data_selected)

  umap_df = as.data.frame(umap_res$layout)
  colnames(umap_df) = c("UMAP1", "UMAP2")
  umap_df$Sample = rownames(umap_df)

  if (color_by == "risk") {
    if (is.null(pred_train) || is.null(pred_test) || is.null(cutoff)) {
      stop("To color by risk, provide pred_train, pred_test, and cutoff.")
    }

    risk_group = rep(NA, nrow(umap_df))
    names(risk_group) = umap_df$Sample

    train_names = rownames(pred_train)
    test_names = rownames(pred_test)

    risk_group[train_names] = ifelse(pred_train > cutoff, "High", "Low")
    risk_group[test_names] = ifelse(pred_test > cutoff, "High", "Low")

    umap_df$risk_group = factor(risk_group[umap_df$Sample], levels = c("Low", "High"))

    # Plot colored by risk stratification
    p = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = risk_group)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_color_manual(values = palette, na.translate = FALSE) +
      labs(title = "UMAP", x = "UMAP1", y = "UMAP2", color = "Risk Group") +
      theme_minimal(base_size = 14)

  } else {
    # Base plot
    p = ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(color = base_color, alpha = 0.7, size = 3) +
      labs(title = "UMAP", x = "UMAP1", y = "UMAP2") +
      theme_minimal(base_size = 14)
  }

  return(p)
}
