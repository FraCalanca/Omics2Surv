#' Bar Chart of Category Counts for a Column
#'
#' Tabulates the values of a single column and draws a labelled bar chart of
#' the counts with \pkg{ggplot2}. Optionally relabels the categories and sets
#' custom fill colours.
#'
#' @param data A data frame containing the column to summarise.
#' @param column Character string naming the column of \code{data} to count.
#' @param labels Optional named character vector remapping the raw category
#'   values (names) to display labels (values). Default is \code{NULL} (raw
#'   values are used as-is).
#' @param fill_colors Optional vector of fill colours passed to
#'   \code{\link[ggplot2]{scale_fill_manual}}. Default is \code{NULL}
#'   (ggplot2 default palette).
#' @param title Optional plot title. Default is \code{NULL}, which uses
#'   \code{paste("Distribution of", column)}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
display_counts <- function(data, column, labels = NULL, fill_colors = NULL, title = NULL) {

  variable <- data[[column]]

  counts <- as.data.frame(table(variable))

  if (!is.null(labels)) {
    counts$status <- factor(counts$variable, levels = names(labels), labels = labels)
  } else {
    counts$status <- as.factor(counts$variable)
  }

  if (is.null(title)) {
    title <- paste("Distribution of", column)
  }

  p <- ggplot2::ggplot(counts, ggplot2::aes(x = status, y = Freq, fill = status)) +
    ggplot2::geom_bar(stat = "identity", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = 3.5, size = 5, color = "black") +
    ggplot2::labs(
      title = title,
      x = column,
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "left")

  if (!is.null(fill_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  }

  return(p)
}


#' Censoring Rate of a Survival Status Vector
#'
#' Computes the proportion of censored observations, i.e. the fraction of
#' entries in \code{status} that are not equal to \code{event_code}.
#'
#' @param status Numeric or integer vector of event indicators.
#' @param event_code Value in \code{status} that denotes an observed event.
#'   Default is \code{1}; every other value is treated as censored.
#'
#' @return A single numeric value in \code{[0, 1]}: the censoring rate.
#'
#' @export
censoring_rate <- function(status, event_code = 1) {
  mean(status != event_code)
}


#' Exploratory PCA and eigencorplot
#'
#' Performs PCA on omics data and generates:
#' 1. 2D PCA plot (PC1 vs PC2)
#' 2. eigencorplot against selected clinical variables
#'
#' @param omics_data Omics data.frame (Linkedomics-like: features × samples)
#' @param clinical_data Clinical data.frame (Linkedomics-like: features × samples)
#' @param metavars Clinical variables for eigencorplot
#'
#' @return list containing:
#' \itemize{
#'   \item pca: PCAtools PCA object
#'   \item biplot: 2D PCA plot
#'   \item eigencorplot: eigencorplot object
#' }
#'
#' @export
explore_pca <- function(
    omics_data,
    clinical_data,
    metavars
) {

  # -----------------------------
  # 1. Clinical data conversion in samples × features
  # -----------------------------
  expr <- omics_data
  clin <- as.data.frame(t(clinical_data))

  # -----------------------------
  # 2. Type conversion
  # -----------------------------
  clin[] <- lapply(clin, function(x) {

    x_num <- suppressWarnings(
      as.numeric(as.character(x))
    )

    if (all(!is.na(x_num))) {
      x_num
    } else {
      as.factor(x)
    }

  })
  # -----------------------------
  # 3. Match samples
  # -----------------------------
  common_ids <- intersect(colnames(expr), rownames(clin))

  expr <- expr[, common_ids, drop = FALSE]
  clin <- clin[common_ids, , drop = FALSE]

  # enforce identical order
  clin <- clin[colnames(expr), , drop = FALSE]

  message(
    "Selected ",
    length(common_ids),
    " matching patients."
  )

  stopifnot(
    identical(colnames(expr), rownames(clin))
  )

  # -----------------------------
  # 4. PCA
  # -----------------------------
  pca_res <- PCAtools::pca(
    expr,
    metadata = clin
  )

  # -----------------------------
  # 5. 2D PCA plot
  # -----------------------------
  biplot <- PCAtools::biplot(
    pca_res,
    x = 1,
    y = 2,
    lab = "",
    singlecol = "",
    legendPosition = "none"
  )

  # -----------------------------
  # 6. Eigencorplot
  # -----------------------------
  eig_plot <- suppressWarnings(
    PCAtools::eigencorplot(
      pca_res,
      metavars = metavars
    )
  )

  return(list(
    pca = pca_res,
    biplot = biplot,
    eigencorplot = eig_plot
  ))
}


#' Collect C-index metrics from a train/test MAE pair
#'
#' Traverses \code{@metadata$analysis} of both MAE objects and assembles a
#' tidy data frame with one row per (model, integration, dataset) combination.
#' Handles all five integration branches (single-omics, early, late, joint,
#' flexynesis) and the extra \code{beta_rho_*} nesting of AFTCoop.
#'
#' @param mae_train MAE object for the training split.
#' @param mae_test  MAE object for the test split.
#' @param aftcoop_rho Numeric scalar selecting which AFTCoop \code{beta_rho}
#'   variant to keep (e.g. \code{0.5}).  Must match one of the keys stored as
#'   \code{beta_rho_<value>} in \code{metrics}.  If \code{NULL} (default) all
#'   three variants are returned as separate rows.
#' @param keep_datasets Character vector of dataset values to retain
#'   (e.g. \code{"mRNA"}, \code{c("mRNA", "early_integrated")},
#'   \code{"multiomics"}).  \code{NULL} (default) keeps all datasets.
#'   Matched against the \code{dataset} column of the output, which mirrors
#'   the keys used in \code{@metadata$analysis}.
#' @param keep_integrations Character vector of integration types to retain
#'   (e.g. \code{"single"}, \code{c("single", "early")}).  \code{NULL}
#'   (default) keeps all.  Matched against the \code{integration} column
#'   (\code{"single"}, \code{"early"}, \code{"late"}, \code{"joint"},
#'   \code{"flexynesis"}).
#'
#' @return A \code{data.frame} with columns \code{model}, \code{integration},
#'   \code{dataset}, \code{harrell_train}, \code{harrell_test}, \code{uno_test},
#'   \code{n_features}.  \code{n_features} is populated for single-omics and
#'   early-integration rows; it is \code{NA} for late, joint and flexynesis
#'   rows.  Other missing values are also \code{NA}.
#'
#' @export
benchmark_models <- function(mae_train, mae_test,
                             aftcoop_rho       = NULL,
                             keep_datasets     = NULL,
                             keep_integrations = NULL) {

  tr <- MultiAssayExperiment::metadata(mae_train)$analysis
  te <- MultiAssayExperiment::metadata(mae_test)$analysis

  results <- list()

  add_row <- function(model, integration, dataset,
                      harrell_train = NA_real_,
                      harrell_test  = NA_real_,
                      uno_test      = NA_real_,
                      n_features    = NA_integer_) {
    results[[length(results) + 1L]] <<- data.frame(
      model = model, integration = integration, dataset = dataset,
      harrell_train = harrell_train, harrell_test = harrell_test,
      uno_test = uno_test, n_features = n_features,
      stringsAsFactors = FALSE
    )
  }

  .nf <- function(x) {
    n <- length(x)
    if (is.character(x) && n > 0L) n else NA_integer_
  }
  .s <- function(x) if (is.numeric(x) && length(x) == 1L) x else NA_real_

  # -----------------------------------------------------------------------
  # 1. STANDARD METHODS (none = single, early = early_*)
  #    tr[[model]][[dataset]]$metrics
  # -----------------------------------------------------------------------
  std_methods <- c("cox_lasso", "cox_adaptive_lasso", "cox_elastic_net", "aft_pen")
  for (model_name in intersect(std_methods, names(tr))) {
    tr_block <- tr[[model_name]]
    te_block <- te[[model_name]]
    for (dataset in names(tr_block)) {
      # early datasets are now tagged "early_<omics>"; everything else is single-omics
      integration <- if (startsWith(dataset, "early_")) "early" else "single"
      add_row(
        model         = model_name,
        integration   = integration,
        dataset       = dataset,
        harrell_train = .s(tr_block[[dataset]]$metrics$harrell_train),
        harrell_test  = .s(te_block[[dataset]]$metrics$harrell_test),
        uno_test      = .s(te_block[[dataset]]$metrics$uno_test),
        n_features    = .nf(tr_block[[dataset]]$train$selected_features)
      )
    }
  }

  # -----------------------------------------------------------------------
  # 2. LATE INTEGRATION
  #    tr$late_integration[[model_method]][[tag]]$metrics
  # -----------------------------------------------------------------------
  if (!is.null(tr$late_integration)) {
    for (late_model in names(tr$late_integration)) {
      tr_lm <- tr$late_integration[[late_model]]
      te_lm <- te$late_integration[[late_model]]
      for (tag in names(tr_lm)) {
        add_row(
          model         = late_model,
          integration   = "late",
          dataset       = tag,
          harrell_train = .s(tr_lm[[tag]]$metrics$harrell_train),
          harrell_test  = .s(te_lm[[tag]]$metrics$harrell_test),
          uno_test      = .s(te_lm[[tag]]$metrics$uno_test)
        )
      }
    }
  }

  # -----------------------------------------------------------------------
  # 3. JOINT MODELS (CoopCox, blockForest, AFTCoop)
  #    tr$joint[[joint_method]][[tag]]$metrics
  #    AFTCoop: metrics has an extra beta_rho_* level
  # -----------------------------------------------------------------------
  if (!is.null(tr$joint)) {
    for (joint_model in names(tr$joint)) {
      tr_jm <- tr$joint[[joint_model]]
      te_jm <- te$joint[[joint_model]]
      for (tag in names(tr_jm)) {
        tr_met <- tr_jm[[tag]]$metrics
        te_met <- te_jm[[tag]]$metrics

        if (joint_model == "AFTCoop") {
          rho_keys <- names(tr_met)
          if (!is.null(aftcoop_rho))
            rho_keys <- rho_keys[rho_keys == paste0("beta_rho_", aftcoop_rho)]
          for (rk in rho_keys) {
            row_model <- if (!is.null(aftcoop_rho)) "AFTCoop" else paste0("AFTCoop_", rk)
            add_row(
              model         = row_model,
              integration   = "joint",
              dataset       = tag,
              harrell_train = .s(tr_met[[rk]]$harrell_train),
              harrell_test  = .s(te_met[[rk]]$harrell_test),
              uno_test      = .s(te_met[[rk]]$uno_test)
            )
          }
        } else {
          add_row(
            model         = joint_model,
            integration   = "joint",
            dataset       = tag,
            harrell_train = .s(tr_met$harrell_train),
            harrell_test  = .s(te_met$harrell_test),
            uno_test      = .s(te_met$uno_test)
          )
        }
      }
    }
  }

  # -----------------------------------------------------------------------
  # 4. FLEXYNESIS
  #    tr$flexynesis[[tag]]$metrics
  # -----------------------------------------------------------------------
  if (!is.null(tr$flexynesis)) {
    for (tag in names(tr$flexynesis)) {
      add_row(
        model         = "flexynesis",
        integration   = "flexynesis",
        dataset       = tag,
        harrell_train = .s(tr$flexynesis[[tag]]$metrics$harrell_train),
        harrell_test  = .s(te$flexynesis[[tag]]$metrics$harrell_test),
        uno_test      = .s(te$flexynesis[[tag]]$metrics$uno_test)
      )
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  out$harrell_train <- round(out$harrell_train, 3)
  out$harrell_test  <- round(out$harrell_test,  3)
  out$uno_test      <- round(out$uno_test,      3)

  if (!is.null(keep_datasets))
    out <- out[out$dataset %in% keep_datasets, , drop = FALSE]
  if (!is.null(keep_integrations))
    out <- out[out$integration %in% keep_integrations, , drop = FALSE]
  if (!nrow(out))
    warning("benchmark_models: all rows were filtered out; check keep_datasets / keep_integrations.")

  out
}


# ---------- defaults (all overridable via c_heatmap args) ----------
STD_METHODS  <- c("cox_lasso","cox_adaptive_lasso","cox_elastic_net","aft_pen")
MODALITY_MAP <- c(mRNA="Omics1", meth="Omics2", early_integrated="Early")
BRANCH_MAP   <- c(late_integration="Late", joint="Joint", flexynesis="Flexynesis")
METRIC_MAP   <- c(harrell_train="H-train", harrell_test="H-test", uno_test="Uno")
SINGLE_OMICS <- c("Omics1","Omics2")
STRAT_ORDER  <- c("Omics1","Omics2","Early","Late","Joint","Flexynesis")
METRIC_ORDER <- c("H-train","H-test","Uno")

# ---------- MAIN ------------------------------------------------------

#' Build a C-index heatmap from a train/test MAE pair
#'
#' Wrapper: harvests the metrics out of a \code{train} and
#' \code{test} MAE object, tidies/classifies them, and returns the faceted
#' heatmap as a \code{ggplot} object.
#'
#' The result is flexible: only the strategies, methods and metrics actually
#' present are drawn. \code{strat_order} and \code{metric_order} are treated as
#' *preferences* — entries not present are skipped, and anything present but
#' not listed is appended (with a message) rather than dropped. So a run with a
#' single omics and a single method (e.g. only \code{cox_lasso} on Omics1)
#' produces a one-row heatmap with no wasted space.
#'
#' @param mae_train,mae_test MAE objects exposing
#'   \code{@metadata[["analysis"]]} for the training and test splits.
#' @param name Character scalar; the dataset label shown as the facet-column
#'   title. Defaults to \code{"Dataset"}.
#' @param aftcoop_rho Numeric scalar; the \code{beta_rho} value to select for
#'   AFTCoop results (e.g. \code{0.3}, \code{0.5}, \code{0.7}).  When
#'   \code{NULL} (default) all available \code{beta_rho_*} sub-keys are kept as
#'   separate rows.  When supplied, only the sub-key whose suffix equals
#'   \code{as.character(aftcoop_rho)} is retained; rows from other \code{rho}
#'   values are silently dropped.
#' @param strat_order Preferred top-to-bottom ordering of strategies. Defaults
#'   to Omics1, Omics2, Early, Late, Joint, Flexynesis. Missing entries are
#'   ignored; unlisted-but-present strategies are appended.
#' @param metric_order Preferred left-to-right ordering of metrics. Defaults to
#'   H-train, H-test, Uno. Handled the same way as \code{strat_order}.
#' @param std_methods,modality_map,branch_map,metric_map,single_omics Lookup
#'   tables used to classify and label paths. Override these if your path keys
#'   or desired labels differ.
#' @param fill_limits Numeric length-2 colour-scale limits. \code{NULL}
#'   (default) derives a padded range from the data.
#' @param text_threshold Value above which tile labels turn white for contrast.
#'   \code{NULL} (default) uses the midpoint of \code{fill_limits}.
#' @param palette ColorBrewer sequential palette name for the fill scale.
#'   Defaults to \code{"YlGnBu"}.
#'
#' @return A \code{ggplot} object. Print it to display, or pass it to
#'   \code{ggsave()} yourself if you want to export it.
#'
#' @export
c_heatmap <- function(mae_train, mae_test, name = "Dataset",
                        aftcoop_rho  = NULL,
                        strat_order  = STRAT_ORDER,
                        metric_order = METRIC_ORDER,
                        std_methods  = STD_METHODS,
                        modality_map = MODALITY_MAP,
                        branch_map   = BRANCH_MAP,
                        metric_map   = METRIC_MAP,
                        single_omics = SINGLE_OMICS,
                        fill_limits  = NULL,
                        text_threshold = NULL,
                        palette      = "YlGnBu") {
  mae_list <- setNames(list(list(train = mae_train, test = mae_test)), name)
  results  <- classify(
    extract_all(mae_list, aftcoop_rho = aftcoop_rho),
    std_methods  = std_methods, modality_map = modality_map,
    branch_map   = branch_map,  metric_map   = metric_map,
    single_omics = single_omics)
  .plot_cindex_heatmap(results, strat_order = strat_order, metric_order = metric_order,
                       single_omics   = single_omics,
                       fill_limits    = fill_limits,
                       text_threshold = text_threshold,
                       palette        = palette)
}

# ---------- INTERNAL: EXTRACTION -------------------------------------

#' Recursively harvest scalar metrics from a nested list
#'
#' Walks an arbitrarily nested list and collects every finite numeric scalar
#' found inside a sub-list named \code{"metrics"}, recording the path taken to
#' reach it.  Handles two depths under \code{metrics}:
#' \itemize{
#'   \item Direct scalars: \code{metrics$harrell_train} etc.
#'   \item One extra level: \code{metrics$beta_rho_0.5$harrell_train} (AFTCoop).
#'     In this case the sub-key (\code{beta_rho_0.5}) is appended to
#'     \code{path} so that \code{classify()} can distinguish the variants.
#' }
#'
#' @param node A list (or any object). Non-list nodes yield an empty result.
#' @param path Character vector of names accumulated on the way down.
#'
#' @return A list of one-row data frames with columns \code{path},
#'   \code{metric} and \code{value}.
#'
#' @keywords internal
#' @noRd
harvest_metrics <- function(node, path = character(0)) {
  out <- list()
  if (!is.list(node)) return(out)
  nms <- names(node); if (is.null(nms)) nms <- as.character(seq_along(node))
  for (i in seq_along(node)) {
    child <- node[[i]]; nm <- nms[i]
    if (identical(nm, "metrics") && is.list(child)) {
      cn <- names(child); if (is.null(cn)) cn <- as.character(seq_along(child))
      for (j in seq_along(child)) {
        v <- child[[j]]
        # normal case: scalar directly under metrics
        if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
          out[[length(out)+1L]] <- data.frame(
            path   = paste(path, collapse = "/"),
            metric = cn[j], value = as.numeric(v),
            stringsAsFactors = FALSE)
          # extra-level case: sub-list (e.g. beta_rho_0.5) with scalars inside
        } else if (is.list(v)) {
          vn <- names(v); if (is.null(vn)) vn <- as.character(seq_along(v))
          for (k in seq_along(v)) {
            vv <- v[[k]]
            if (is.numeric(vv) && length(vv) == 1L && is.finite(vv))
              out[[length(out)+1L]] <- data.frame(
                path   = paste(c(path, cn[j]), collapse = "/"),
                metric = vn[k], value = as.numeric(vv),
                stringsAsFactors = FALSE)
          }
        }
      }
    } else out <- c(out, harvest_metrics(child, c(path, nm)))
  }
  out
}

#' Extract metrics from a single MAE object, optionally filtering AFTCoop rho
#'
#' Runs \code{harvest_metrics()} on the \code{analysis} entry of an MAE
#' object's metadata, optionally retains only the selected \code{beta_rho_*}
#' sub-key for AFTCoop rows, and tags the result with a dataset label.
#'
#' @param mae An MAE object exposing \code{mae@metadata[["analysis"]]}.
#' @param ds Character scalar; the dataset name to attach.
#' @param aftcoop_rho Numeric or \code{NULL}; see \code{c_heatmap}.
#'
#' @return A data frame with columns \code{path}, \code{metric}, \code{value},
#'   \code{dataset}, or \code{NULL} if no metrics were found after filtering.
#'
#' @keywords internal
#' @noRd
extract_one <- function(mae, ds, aftcoop_rho = NULL) {
  r <- harvest_metrics(mae@metadata[["analysis"]])
  if (!length(r)) return(NULL)
  d <- do.call(rbind, r)
  d$dataset <- ds
  # filter AFTCoop rho if requested
  if (!is.null(aftcoop_rho)) {
    rho_tag   <- paste0("beta_rho_", aftcoop_rho)   # e.g. "beta_rho_0.5"
    # rows with a beta_rho_* suffix in their path: keep only the chosen one
    has_rho   <- grepl("beta_rho_", d$path, fixed = TRUE)
    keep_rho  <- grepl(rho_tag,     d$path, fixed = TRUE)
    d <- d[!has_rho | keep_rho, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    # strip the beta_rho_* segment from the path so classify() sees a clean path
    d$path <- sub(paste0("/", rho_tag), "", d$path, fixed = TRUE)
  }
  d
}

#' Extract metrics from a list of train/test MAE pairs
#'
#' @param mae_list A named list where each element is a list with \code{train}
#'   and \code{test} MAE objects.
#' @param aftcoop_rho Passed through to \code{extract_one}.
#'
#' @return A single data frame combining metrics across all datasets/splits.
#'
#' @keywords internal
#' @noRd
extract_all <- function(mae_list, aftcoop_rho = NULL) {
  res <- list()
  for (ds in names(mae_list)) {
    res[[length(res)+1L]] <- extract_one(mae_list[[ds]]$train, ds, aftcoop_rho)
    res[[length(res)+1L]] <- extract_one(mae_list[[ds]]$test,  ds, aftcoop_rho)
  }
  do.call(rbind, res)
}

#' Annotate harvested metrics with strategy, method, omics and tidy labels
#'
#' Parses the slash-separated \code{path} of each metric to decide whether it
#' came from a standard single-/early-integration method (\code{std_methods})
#' or from a branch, then derives \code{strategy}, \code{method}, \code{omics}
#' and tidy \code{metric} labels from the supplied lookup tables. Paths it can't
#' classify get an \code{NA} strategy and are later dropped by the plotter.
#'
#' @param df A data frame from \code{extract_all()} with \code{path},
#'   \code{metric}, \code{value} columns.
#' @param std_methods,modality_map,branch_map,metric_map,single_omics Lookup
#'   tables (see \code{c_heatmap}).
#'
#' @return A data frame with columns \code{dataset}, \code{strategy},
#'   \code{method}, \code{omics}, \code{metric}, \code{value}, \code{path}.
#'
#' @keywords internal
#' @noRd
classify <- function(df,
                     std_methods  = STD_METHODS,
                     modality_map = MODALITY_MAP,
                     branch_map   = BRANCH_MAP,
                     metric_map   = METRIC_MAP,
                     single_omics = SINGLE_OMICS) {
  parts  <- strsplit(df$path, "/", fixed = TRUE)
  first  <- vapply(parts, `[`, "", 1)
  second <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, "")
  is_std   <- first %in% std_methods
  # Early integration is now stored under a per-run tag ("early_<omics>")
  # rather than a fixed "early_integrated" key; detect it by prefix so the
  # rows are labelled "Early" regardless of which omics were fused.
  is_early <- is_std & !is.na(second) & startsWith(second, "early_")
  df$strategy <- ifelse(is_early, "Early",
                        ifelse(is_std, unname(modality_map[second]),
                               unname(branch_map[first])))
  df$method   <- ifelse(is_std, first,
                        vapply(parts, function(p) paste(p[-1], collapse = "/"), ""))
  df$metric   <- ifelse(df$metric %in% names(metric_map),
                        unname(metric_map[df$metric]), df$metric)
  df$omics    <- ifelse(df$strategy %in% single_omics, "Single-omics", "Multi-omics")
  df[, c("dataset","strategy","method","omics","metric","value","path")]
}

# ---------- INTERNAL: HEATMAP ----------------------------------------

#' Keep the requested ordering for present levels and append the rest sorted
#' @keywords internal
#' @noRd
.order_levels <- function(present, preferred, what = "level") {
  present <- unique(present)
  extra   <- setdiff(present, preferred)
  if (length(extra))
    message("plot_heatmap: ", what, "(s) not in the requested order, appended: ",
            paste(extra, collapse = ", "))
  c(intersect(preferred, present), sort(extra))
}

#' Build a faceted C-index heatmap (adapts to the data present)
#'
#' \code{geom_tile} heatmap with datasets in facet columns, models in rows and
#' metrics on the x-axis.
#'
#' @param results A classified data frame from \code{classify()}.
#' @param strat_order,metric_order Preferred orderings (see \code{c_heatmap}).
#' @param single_omics Strategy labels counted as single-omics (for the divider).
#' @param fill_limits Length-2 colour limits, or \code{NULL} to derive from data.
#' @param text_threshold Label-whitening cutoff, or \code{NULL} for the midpoint.
#' @param palette ColorBrewer sequential palette name (default \code{"YlGnBu"}).
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
#' @noRd
.plot_cindex_heatmap <- function(results,
                                 strat_order    = STRAT_ORDER,
                                 metric_order   = METRIC_ORDER,
                                 single_omics   = SINGLE_OMICS,
                                 fill_limits    = NULL,
                                 text_threshold = NULL,
                                 palette        = "YlGnBu") {
  results <- results |> dplyr::filter(!is.na(strategy), !is.na(value))
  if (!nrow(results)) stop("plot_heatmap: no classifiable metrics to plot.")

  results$strategy <- factor(results$strategy,
                             levels = .order_levels(results$strategy, strat_order, "strategy"))
  results$metric   <- factor(results$metric,
                             levels = .order_levels(results$metric, metric_order, "metric"))
  results$dataset  <- factor(results$dataset)

  row_tab <- results |> dplyr::distinct(strategy, method, omics) |>
    dplyr::arrange(strategy, method) |> dplyr::mutate(model = paste(strategy, method))
  results <- results |> dplyr::mutate(model = paste(strategy, method))
  results$model <- factor(results$model, levels = rev(row_tab$model))

  if (is.null(fill_limits)) {
    rng <- range(results$value); pad <- max(0.01, diff(rng) * 0.05)
    fill_limits <- c(floor((rng[1]-pad)*100)/100, ceiling((rng[2]+pad)*100)/100)
  }
  if (is.null(text_threshold)) text_threshold <- mean(fill_limits)

  p <- ggplot2::ggplot(results, ggplot2::aes(metric, model, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value), colour = value > text_threshold),
                       size = 2.7, show.legend = FALSE) +
    ggplot2::facet_grid(. ~ dataset) +
    ggplot2::scale_fill_distiller(palette = palette, direction = 1,
                                  limits = fill_limits, name = "C-index") +
    ggplot2::scale_colour_manual(values = c(`TRUE`="white", `FALSE`="grey15")) +
    ggplot2::labs(title = "C-index across models, datasets and metrics",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", size = 13),
                   axis.text.x = ggplot2::element_text(size = 8))

  n_multi  <- sum(row_tab$omics == "Multi-omics")
  n_single <- sum(row_tab$omics == "Single-omics")
  if (n_multi > 0 && n_single > 0)
    p <- p + ggplot2::geom_hline(yintercept = n_multi + 0.5, linewidth = 1.1, colour = "black")

  p
}


#' @keywords internal
#' @noRd
.risk_colors <- function() {
  c("High risk" = "red", "Low risk" = "deepskyblue")
}

#' @keywords internal
#' @noRd
.check_strata <- function(x) {
  if (!is.data.frame(x)) {
    stop("Expected a stratification data frame.", call. = FALSE)
  }
  miss <- setdiff(c("score", "time", "status", "risk_group"), names(x))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "),
         ". Use patient_stratification().", call. = FALSE)
  }
  invisible(TRUE)
}

# Accept either the patient_stratification() list or a plain strata data frame.
#' @keywords internal
#' @noRd
.as_strata <- function(x) {
  if (is.list(x) && !is.data.frame(x) && !is.null(x$strata)) x$strata else x
}

# Locate a model node in metadata()$analysis. With `container` given, look at
# analysis[[container]][[model]] (needed to disambiguate e.g. late_integration,
# whose sub-models share names with the top-level ones). Otherwise match at the
# top level, then one level down into container entries.
#' @keywords internal
#' @noRd
.locate_model <- function(an, model, container = NULL) {
  if (!is.null(container)) {
    if (!container %in% names(an)) {
      stop(sprintf("Container '%s' not found in metadata(mae)$analysis.",
                   container), call. = FALSE)
    }
    sub <- an[[container]]
    if (!model %in% names(sub)) {
      stop(sprintf("Model '%s' not found in container '%s'. Available: %s",
                   model, container, paste(names(sub), collapse = ", ")),
           call. = FALSE)
    }
    return(sub[[model]])
  }
  if (model %in% names(an)) return(an[[model]])
  for (cont in names(an)) {
    sub <- an[[cont]]
    if (is.list(sub) && model %in% names(sub)) return(sub[[model]])
  }
  stop(sprintf("Model '%s' not found in metadata(mae)$analysis.", model),
       call. = FALSE)
}

# Resolve which split ("train"/"test") to use in a node. Only "train"/"test"
# count as splits; siblings such as "metrics" or "late_method" are ignored.
#' @keywords internal
#' @noRd
.resolve_split <- function(node, split = "auto") {
  splits <- intersect(c("train", "test"), names(node))
  if (identical(split, "auto")) {
    if (length(splits) != 1) {
      stop("Ambiguous or missing split; set `split` to one of: ",
           paste(c(splits, "train", "test"), collapse = ", "), call. = FALSE)
    }
    return(splits)
  }
  if (!split %in% names(node)) {
    stop(sprintf("Split '%s' not found. Available: %s",
                 split, paste(splits, collapse = ", ")), call. = FALSE)
  }
  split
}

# From a located model node, return the container that holds the splits
# (train/test/metrics). Cox models nest results under a modality level; joint
# sub-models (blockForest, CoopCox, AFTCoop) and flexynesis do not.
#' @keywords internal
#' @noRd
.model_split_node <- function(node, modality = NULL) {
  if (any(c("train", "test") %in% names(node))) {
    return(node)                                   # no modality level
  }
  keys <- setdiff(names(node), "metrics")
  if (is.null(modality)) {
    # A single per-modality entry (e.g. one "early_<omics>" or "late_<omics>"
    # tag) is unambiguous, so resolve it automatically; otherwise ask.
    if (length(keys) == 1L) return(node[[keys]])
    stop("This model has per-modality results; supply `modality` (one of: ",
         paste(keys, collapse = ", "), ").",
         call. = FALSE)
  }
  if (!modality %in% names(node)) {
    stop(sprintf("Modality '%s' not found. Available: %s",
                 modality, paste(keys, collapse = ", ")),
         call. = FALSE)
  }
  node[[modality]]
}

# The assay that holds a model's features: the modality itself when given, or
# the integrated assay for joint models (modality NULL). The integrated assay
# is now named with a per-run tag ("early_<omics>"), so when no modality is
# supplied we resolve the single "early_*" assay present. Shared by plot_umap,
# plot_heatmap and enrich_features so the feature space, projection space and
# enrichment background are resolved the same way.
#' @keywords internal
#' @noRd
.default_experiment <- function(mae, modality = NULL) {
  if (!is.null(modality)) return(modality)
  ex    <- names(mae)
  early <- grep("^early_", ex, value = TRUE)
  if (length(early) == 1L) return(early)
  if ("early_integrated" %in% ex) return("early_integrated")   # legacy key
  if (length(early) > 1L) {
    stop("Multiple early-integration assays found (",
         paste(early, collapse = ", "), "); specify `experiment`.",
         call. = FALSE)
  }
  stop("Cannot infer the experiment to draw; specify `experiment` or ",
       "`modality`.", call. = FALSE)
}

# Find the variable-importance object within a model node.
#' @keywords internal
#' @noRd
.get_varimp <- function(node) {
  # Importance may sit at this level, one level down (the per-run tag added by
  # the joint / flexynesis storage), or two levels down (the train/test split),
  # either directly ($variable_importance / $feature_importance) or under
  # $extra (blockForest). Collect the first match found at any of those depths.
  found <- NULL
  .pick <- function(x) {
    if (!is.list(x) || !is.null(found)) return(invisible())
    if (!is.null(x$variable_importance)) { found <<- x$variable_importance; return(invisible()) }
    if (!is.null(x$feature_importance))  { found <<- x$feature_importance;  return(invisible()) }
    ex <- x$extra
    if (!is.null(ex$variable_importance)) { found <<- ex$variable_importance; return(invisible()) }
    if (!is.null(ex$feature_importance))  { found <<- ex$feature_importance }
  }

  .pick(node)                                        # model level
  if (is.null(found)) {
    for (s in setdiff(names(node), "metrics")) {     # tag level
      child <- node[[s]]
      .pick(child)
      if (is.null(found) && is.list(child)) {        # train / test split
        for (s2 in setdiff(names(child), "metrics")) .pick(child[[s2]])
      }
      if (!is.null(found)) break
    }
  }

  if (!is.null(found)) return(found)
  stop("No variable importance found for this model. `top_features()` is ",
       "for models without built-in feature selection (e.g. blockForest, ",
       "flexynesis).", call. = FALSE)
}

# Coerce a variable-importance object (named vector or data frame) to a named
# numeric vector (names = feature, value = importance). For flexynesis-style
# tables the name is prefixed with its layer (e.g. "meth_AGPAT4") so it matches
# an integrated assay and layers are kept distinct. Duplicated keys are
# collapsed by max.
#' @keywords internal
#' @noRd
.varimp_to_named <- function(x) {
  if (is.data.frame(x)) {
    if (all(c("name", "importance") %in% names(x))) {           # flexynesis
      key <- if ("layer" %in% names(x)) paste(x$layer, x$name, sep = "_")
      else x$name
      return(tapply(x$importance, key, max, na.rm = TRUE))
    }
    v <- as.numeric(x[[2]]); names(v) <- as.character(x[[1]])   # generic
    return(tapply(v, names(v), max, na.rm = TRUE))
  }
  if (is.numeric(x) && !is.null(names(x))) return(x)
  stop("Unrecognised variable-importance format.", call. = FALSE)
}

# Resolve a model's feature set: selected_features() for penalised models,
# top_features() for importance-based ones (blockForest, flexynesis).
#' @keywords internal
#' @noRd
.resolve_model_features <- function(src, model, modality, top_feat, features,
                                    rho = NULL) {
  if (!is.null(features)) return(features)
  if (model %in% c("blockForest", "flexynesis")) {
    if (is.null(top_feat)) top_features(src, model)
    else                   top_features(src, model, n = top_feat)
  } else {
    selected_features(src, model, modality, rho = rho)
  }
}

# Feature x sample matrix restricted to the given features and to the
# stratified samples, with the risk-group factor aligned to the columns.
#' @keywords internal
#' @noRd
.feature_matrix <- function(mae, strat, features, experiment, assay) {
  se  <- mae[[experiment]]
  mat <- SummarizedExperiment::assay(se, assay)   # features x samples

  present <- intersect(features, rownames(mat))
  if (!length(present)) {
    stop("None of the model features are present in experiment '",
         experiment, "'.", call. = FALSE)
  }
  n_missing <- length(unique(features)) - length(present)
  if (n_missing > 0) {
    warning(n_missing, " model feature(s) not found in the assay and skipped.",
            call. = FALSE)
  }
  samp <- intersect(colnames(mat), strat$sample_id)
  if (!length(samp)) {
    stop("No stratified samples found in experiment '", experiment, "'.",
         call. = FALSE)
  }
  list(mat    = mat[present, samp, drop = FALSE],
       groups = droplevels(strat$risk_group[match(samp, strat$sample_id)]))
}

# Extract predictions + survival from a MAE as a tidy data frame.
#' @keywords internal
#' @noRd
.extract_predictions <- function(mae, model, modality = NULL, split = "auto",
                                 container = NULL,
                                 rho = NULL,
                                 time_col = "overall_survival",
                                 status_col = "status") {
  if (!methods::is(mae, "MultiAssayExperiment")) {
    stop("`mae` must be a MultiAssayExperiment.", call. = FALSE)
  }
  an <- S4Vectors::metadata(mae)$analysis
  if (is.null(an)) stop("metadata(mae)$analysis not found.", call. = FALSE)

  node  <- .model_split_node(.locate_model(an, model, container), modality)
  split <- .resolve_split(node, split)

  snode <- node[[split]]
  pred  <- snode$pred
  if (is.null(pred)) pred <- snode

  # AFTCoop: pred is a named list (one vector per rho); resolve to one vector.
  pred <- .resolve_aftcoop_rho(pred, rho, what = "pred")

  if (is.matrix(pred) || is.array(pred)) {
    ids <- rownames(pred); score <- as.numeric(pred[, 1])
  } else if (is.data.frame(pred)) {
    ids <- rownames(pred); score <- as.numeric(pred[[1]])
  } else {
    ids <- names(pred);    score <- as.numeric(pred)
  }

  # Late integration stores an unnamed consensus rank alongside a merged_risks
  # table carrying patient_id; use it as the sample identifiers.
  if (is.null(ids) && is.data.frame(snode$merged_risks) &&
      "patient_id" %in% names(snode$merged_risks) &&
      nrow(snode$merged_risks) == length(score)) {
    ids <- as.character(snode$merged_risks$patient_id)
  }

  cd <- as.data.frame(MultiAssayExperiment::colData(mae),
                      stringsAsFactors = FALSE)

  if (is.null(ids)) {
    stop("Predictions for '", model, "' have no sample IDs; cannot align to ",
         "colData. Store the predictions with sample barcodes as names.",
         call. = FALSE)
  }

  idx <- match(ids, rownames(cd))
  df  <- data.frame(
    sample_id = as.character(ids),
    score     = score,
    time      = as.numeric(cd[[time_col]][idx]),
    status    = as.integer(cd[[status_col]][idx]),
    stringsAsFactors = FALSE
  )
  df <- df[stats::complete.cases(df[, c("score", "time", "status")]), ,
           drop = FALSE]
  rownames(df) <- NULL
  df
}

# Compute a cut-point from the training scores.
#' @keywords internal
#' @noRd
.find_threshold <- function(df, method, threshold, minprop) {
  switch(
    method,
    median = stats::median(df$score, na.rm = TRUE),
    user = {
      if (is.null(threshold) || !is.numeric(threshold)) {
        stop("`threshold` must be numeric when method = 'user'.",
             call. = FALSE)
      }
      threshold[1]
    },
    survcutpoint = {
      if (!requireNamespace("survminer", quietly = TRUE)) {
        stop("Package 'survminer' is required for method = 'survcutpoint'.",
             call. = FALSE)
      }
      cp <- survminer::surv_cutpoint(df, time = "time", event = "status",
                                     variables = "score", minprop = minprop)
      as.numeric(summary(cp)$cutpoint)
    }
  )
}

# Assign risk groups given a threshold (higher score = higher risk).
#' @keywords internal
#' @noRd
.stratify <- function(df, threshold, inclusive) {
  hi <- if (inclusive) df$score >= threshold else df$score > threshold
  df$risk_group <- factor(ifelse(hi, "High risk", "Low risk"),
                          levels = c("High risk", "Low risk"))
  attr(df, "threshold") <- as.numeric(threshold)
  df
}

# Histogram of the training scores, coloured by group, with the threshold.
#' @keywords internal
#' @noRd
.score_histogram <- function(df, threshold, colors) {
  bins <- grDevices::nclass.Sturges(df$score)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$score)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = .data$risk_group),
      bins = bins, colour = "white", alpha = 0.9, position = "identity") +
    ggplot2::scale_fill_manual(values = colors, name = NULL) +
    ggplot2::geom_vline(xintercept = threshold,
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(x = "Risk score (training)", y = "Count") +
    ggplot2::theme_bw()
}

# ---------------------------------------------------------------------------
# 1. patient_stratification
# ---------------------------------------------------------------------------

#' Stratify Patients into Risk Groups from Training and Test MAEs
#'
#' Single entry point that (1) extracts a model's per-patient risk scores from
#' the training \code{MultiAssayExperiment}, (2) derives a cut-point on those
#' training scores, and (3) extracts and stratifies the test predictions using
#' that cut-point. Survival outcome (\code{overall_survival} and \code{status})
#' is read directly from \code{colData(mae)}; extraction and coercion are handled
#' internally. Predictions are read from
#' \code{metadata(mae)$analysis[[model]][[modality]][[split]]$pred}.
#'
#' @param mae_train A \code{MultiAssayExperiment} holding the training split;
#'   the threshold and the optional score histogram are derived from it.
#'   \code{colData} must contain \code{"overall_survival"} and \code{"status"}.
#' @param mae_test A \code{MultiAssayExperiment} holding the test split, to be
#'   stratified with the training threshold. If \code{NULL}, the threshold is
#'   applied back to the training set.
#' @param model Character string naming the model in
#'   \code{metadata()$analysis} (e.g. \code{"cox_elastic_net"}).
#' @param modality Character string naming the omics modality
#'   (e.g. an \code{"early_<omics>"} tag, \code{"mRNA"}, \code{"meth"}). Leave
#'   \code{NULL} for models stored without a modality level (joint methods such
#'   as \code{"blockForest"}, and \code{"flexynesis"}), or to auto-resolve a
#'   model that has a single per-modality entry (e.g. one early-integration
#'   tag, or the single ensemble under \code{container = "late_integration"}).
#' @param container Optional container to disambiguate models whose names are
#'   reused across integration strategies. Use \code{"late_integration"} to
#'   read the consensus prediction of a late-integration ensemble (e.g.
#'   \code{model = "cox_lasso"}, \code{container = "late_integration"}); the
#'   sample IDs are taken from its \code{merged_risks} table.
#' @param method Threshold method:
#'   \itemize{
#'     \item \code{"median"} (alias \code{"mediana"}): median of the training
#'       scores.
#'     \item \code{"user"}: user-supplied value, passed via \code{threshold}.
#'     \item \code{"survcutpoint"}: optimal cut-point via
#'       \code{\link[survminer]{surv_cutpoint}}.
#'   }
#' @param rho Numeric. For \code{model = "AFTCoop"} only: selects the
#'   \code{beta_rho_<rho>} variant (e.g. \code{rho = 0.5}). If \code{NULL}
#'   (default) and only one rho was trained, it is resolved automatically;
#'   if multiple rho values are present, an error lists the available keys.
#' @param threshold Numeric cut-point. Required when \code{method = "user"}.
#' @param plot_scores Logical. If \code{TRUE}, draws and returns a histogram of
#'   the training scores coloured by risk group, with the threshold marked.
#'   Default is \code{FALSE}.
#' @param minprop Minimum group proportion for
#'   \code{\link[survminer]{surv_cutpoint}}. Default is 0.1.
#' @param inclusive Logical. If \code{TRUE}, scores equal to the threshold are
#'   assigned to \code{"High risk"}. Default is \code{FALSE}.
#' @param colors Named character vector mapping \code{"High risk"} and
#'   \code{"Low risk"} to colours. Default is
#'   \code{c("High risk" = "red", "Low risk" = "deepskyblue")}.
#' @param time_col,status_col Names of the survival time and event columns in
#'   \code{colData}. Default to \code{"overall_survival"} and \code{"status"}.
#'
#' @return A \code{list} with:
#'   \itemize{
#'     \item \code{strata}: the stratified test \code{data.frame} (pass to
#'       \code{\link{plot_km}} / \code{\link{plot_umap}}).
#'     \item \code{threshold}: the numeric cut-point used.
#'     \item \code{method}: the threshold method applied.
#'     \item \code{train}: the stratified training \code{data.frame}.
#'     \item \code{histogram}: a \code{ggplot} object, or \code{NULL}.
#'   }
#'
#' @seealso \code{\link{plot_km}}, \code{\link{plot_umap}}
#'
#' @export
patient_stratification <- function(mae_train,
                                   mae_test = NULL,
                                   model,
                                   modality = NULL,
                                   container = NULL,
                                   rho = NULL,
                                   method = c("median", "mediana",
                                              "user", "survcutpoint"),
                                   threshold = NULL,
                                   plot_scores = FALSE,
                                   minprop = 0.1,
                                   inclusive = FALSE,
                                   colors = .risk_colors(),
                                   time_col = "overall_survival",
                                   status_col = "status") {
  method <- match.arg(method)
  if (method == "mediana") method <- "median"

  # AFT models store a linear predictor oriented the opposite way (higher =
  # longer survival), so negate their score to enforce "higher = higher risk".
  invert <- .is_aft_model(model)

  train_df <- .extract_predictions(mae_train, model, modality,
                                   container = container, rho = rho,
                                   time_col = time_col, status_col = status_col)
  if (invert) train_df$score <- -train_df$score

  thr      <- .find_threshold(train_df, method, threshold, minprop)
  train_st <- .stratify(train_df, thr, inclusive)

  if (is.null(mae_test)) {
    test_st <- train_st
  } else {
    test_df <- .extract_predictions(mae_test, model, modality,
                                    container = container, rho = rho,
                                    time_col = time_col, status_col = status_col)
    if (invert) test_df$score <- -test_df$score
    test_st <- .stratify(test_df, thr, inclusive)
  }

  hist_plot <- NULL
  if (isTRUE(plot_scores)) {
    hist_plot <- .score_histogram(train_st, thr, colors)
    print(hist_plot)
  }

  list(strata = test_st, threshold = as.numeric(thr), method = method,
       train = train_st, histogram = hist_plot)
}

# ---------------------------------------------------------------------------
# 2. plot_km
# ---------------------------------------------------------------------------

#' Kaplan-Meier Plot by Risk Group
#'
#' Draws Kaplan-Meier survival curves for the \code{"High risk"} and
#' \code{"Low risk"} groups produced by \code{\link{patient_stratification}},
#' using \code{\link[survminer]{ggsurvplot}}.
#'
#' @param x Either the \code{list} returned by
#'   \code{\link{patient_stratification}} or its \code{strata}
#'   \code{data.frame}.
#' @param colors Named character vector mapping risk groups to colours. Default
#'   is \code{c("High risk" = "red", "Low risk" = "deepskyblue")}.
#' @param pval Logical, passed to \code{\link[survminer]{ggsurvplot}}. Default
#'   is \code{TRUE}.
#' @param conf_int Logical, passed to \code{\link[survminer]{ggsurvplot}} as
#'   \code{conf.int} to draw the confidence bands. Default is \code{TRUE}.
#' @param risk_table Logical, passed to \code{\link[survminer]{ggsurvplot}} as
#'   \code{risk.table}. Default is \code{TRUE}.
#' @param time_unit Character string used as the x-axis label. Default is
#'   \code{"Time"}.
#' @param ... Additional arguments passed to
#'   \code{\link[survminer]{ggsurvplot}}.
#'
#' @return A \code{ggsurvplot} object.
#'
#' @seealso \code{\link{patient_stratification}}, \code{\link{plot_umap}}
#'
#' @export
plot_km <- function(x,
                    colors = .risk_colors(),
                    pval = TRUE,
                    conf_int = TRUE,
                    risk_table = TRUE,
                    time_unit = "Time",
                    ...) {
  strat <- .as_strata(x)
  .check_strata(strat)
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required for plot_km().", call. = FALSE)
  }
  strat <- droplevels(strat)
  fit <- survival::survfit(
    survival::Surv(time, status) ~ risk_group, data = strat)
  survminer::ggsurvplot(
    fit, data = strat,
    palette      = unname(colors[levels(strat$risk_group)]),
    pval         = pval,
    conf.int     = conf_int,
    risk.table   = risk_table,
    legend.title = "",
    legend.labs  = levels(strat$risk_group),
    xlab         = time_unit,
    ...
  )
}

# ---------------------------------------------------------------------------
# 3. top_features
# ---------------------------------------------------------------------------

#' Top Features by Variable Importance
#'
#' Extracts a model's variable importance from
#' \code{metadata(mae)$analysis}, ranks it in decreasing order and cuts it to a
#' feature set. Intended for models that do \emph{not} perform built-in feature
#' selection (e.g. \code{"blockForest"}, \code{"flexynesis"}); the returned
#' features can then be passed to \code{\link{plot_umap}}. The model node is
#' located automatically, including methods nested under container entries such
#' as \code{"joint"} (e.g. \code{blockForest} at
#' \code{analysis$joint$blockForest}).
#'
#' Both storage formats are handled: a named numeric importance vector, and a
#' \code{data.frame} (either the \code{flexynesis} layout with \code{name} and
#' \code{importance} columns, or a two-column feature/importance table).
#' Duplicated feature names (e.g. a feature reported across several
#' \code{flexynesis} layers) are collapsed by their maximum importance.
#'
#' @param mae A \code{MultiAssayExperiment} whose
#'   \code{metadata()$analysis} holds the fitted model.
#' @param model Character string naming the model (e.g. \code{"blockForest"},
#'   \code{"flexynesis"}).
#' @param n Integer. Number of top features to keep. Default is 50. Set to
#'   \code{NULL} to keep all features passing \code{min_importance}.
#' @param min_importance Numeric. If supplied, only features with importance
#'   greater than or equal to this value are kept (applied before \code{n}).
#'   Default is \code{NULL}.
#'
#' @return A character vector of feature names, ordered by decreasing
#'   importance. Importance values are attached as the \code{"importance"}
#'   attribute.
#'
#' @seealso \code{\link{plot_umap}}
#'
#' @export
top_features <- function(mae, model, n = 50, min_importance = NULL) {
  if (!methods::is(mae, "MultiAssayExperiment")) {
    stop("`mae` must be a MultiAssayExperiment.", call. = FALSE)
  }
  an <- S4Vectors::metadata(mae)$analysis
  if (is.null(an)) stop("metadata(mae)$analysis not found.", call. = FALSE)

  node <- .locate_model(an, model)
  imp  <- sort(.varimp_to_named(.get_varimp(node)), decreasing = TRUE)

  if (!is.null(min_importance)) imp <- imp[imp >= min_importance]
  if (!is.null(n))              imp <- utils::head(imp, n)

  feats <- names(imp)
  attr(feats, "importance") <- unname(imp)
  feats
}

# ---------------------------------------------------------------------------
# 4. selected_features
# ---------------------------------------------------------------------------

#' Selected Features of a Penalised Model
#'
#' Returns the features retained by a model that performs built-in feature
#' selection (e.g. the non-zero-coefficient predictors of \code{"cox_lasso"},
#' \code{"cox_adaptive_lasso"}, \code{"cox_elastic_net"}), read from
#' \code{metadata(mae)$analysis[[model]][[modality]][[split]]$selected_features}.
#' The returned features can be passed to \code{\link{plot_umap}}. For models
#' \emph{without} built-in selection (\code{"blockForest"},
#' \code{"flexynesis"}), use \code{\link{top_features}} instead.
#'
#' @param mae A \code{MultiAssayExperiment} whose \code{metadata()$analysis}
#'   holds the fitted model.
#' @param model Character string naming the model (e.g. \code{"cox_lasso"}).
#' @param modality Character string naming the omics modality
#'   (e.g. an \code{"early_<omics>"} tag, \code{"mRNA"}, \code{"meth"}). Leave
#'   \code{NULL} for models stored without a modality level (joint methods), or
#'   to auto-resolve a model that has a single per-modality entry.
#' @param split Which split to read: \code{"auto"} (default; the single split
#'   present besides \code{"metrics"}, i.e. \code{"train"} in a training MAE)
#'   or an explicit \code{"train"} / \code{"test"}.
#' @param rho Numeric. For \code{model = "AFTCoop"} only: selects the
#'   \code{beta_rho_<rho>} feature set (e.g. \code{rho = 0.5}). If \code{NULL}
#'   (default) and only one rho was trained, it is resolved automatically.
#'
#' @return A character vector of selected feature names.
#'
#' @seealso \code{\link{top_features}}, \code{\link{plot_umap}}
#'
#' @export
selected_features <- function(mae, model, modality = NULL, split = "auto",
                              rho = NULL) {
  if (!methods::is(mae, "MultiAssayExperiment")) {
    stop("`mae` must be a MultiAssayExperiment.", call. = FALSE)
  }
  an    <- S4Vectors::metadata(mae)$analysis
  if (is.null(an)) stop("metadata(mae)$analysis not found.", call. = FALSE)
  node  <- .model_split_node(.locate_model(an, model), modality)
  split <- .resolve_split(node, split)

  sf <- node[[split]]$selected_features
  if (is.null(sf)) {
    stop("No `selected_features` for model '", model,
         "'. Use top_features() for importance-based models ",
         "(blockForest, flexynesis).", call. = FALSE)
  }
  # AFTCoop: selected_features is a list keyed by rho; resolve to one vector.
  sf <- .resolve_aftcoop_rho(sf, rho, what = "selected_features")
  if (is.data.frame(sf)) {
    if ("name" %in% names(sf)) {
      sf <- if ("layer" %in% names(sf)) paste(sf$layer, sf$name, sep = "_")
      else sf$name
    } else {
      col <- intersect("feature", names(sf))
      sf  <- if (length(col)) sf[[col[1]]] else sf[[1]]
    }
  } else if (is.list(sf)) {
    sf <- unlist(sf, use.names = FALSE)
  }
  unique(as.character(sf))
}

# ---------------------------------------------------------------------------
# 5. plot_umap
# ---------------------------------------------------------------------------

#' UMAP of Patients Coloured by Risk Group
#'
#' Projects patients into two dimensions with \code{\link[uwot]{umap}} and
#' colours them by the risk group from \code{\link{patient_stratification}}.
#' The projection is restricted to the features of the given \code{model}:
#' \code{selected_features()} for models with built-in selection (Cox penalised
#' models) and \code{top_features()} for importance-based models
#' (\code{"blockForest"}, \code{"flexynesis"}); the choice is made from
#' \code{model}. Assay samples are matched to stratified patients by barcode.
#'
#' Feature information (selected features / variable importance) is a training
#' artifact, so it is read from \code{mae_train}. The projection itself is drawn
#' on \code{mae} (typically the stratified test set).
#'
#' @param mae A \code{MultiAssayExperiment} providing the feature matrix to
#'   project (typically the stratified test object).
#' @param x Either the \code{list} returned by
#'   \code{\link{patient_stratification}} or its \code{strata}
#'   \code{data.frame}.
#' @param model Character string naming the model whose features drive the
#'   projection (e.g. \code{"cox_lasso"}, \code{"blockForest"},
#'   \code{"flexynesis"}).
#' @param modality Character string naming the omics modality. It locates a
#'   penalised model's selected features (ignored for \code{"blockForest"} /
#'   \code{"flexynesis"}) and, unless \code{experiment} is given, also selects
#'   the assay to project on. Default is \code{NULL}: the single per-run
#'   integration tag (\code{"early_<omics>"}) is then resolved automatically;
#'   pass the tag (or an omic name) explicitly to disambiguate.
#' @param top_feat Integer number of top features to keep. Used \emph{only} when
#'   \code{model} is \code{"blockForest"} or \code{"flexynesis"}; ignored
#'   otherwise. If \code{NULL} (default), \code{\link{top_features}}'s own
#'   default is used.
#' @param mae_train A \code{MultiAssayExperiment} from which to read the model's
#'   features. Defaults to \code{mae}; set this when feature information lives in
#'   the training object.
#' @param features Optional character vector of feature (row) names. If
#'   supplied, it overrides the automatic retrieval from \code{model}.
#' @param experiment Name or index of the experiment to project on. If
#'   \code{NULL} (default), the assay matching \code{modality} is used, or the
#'   single \code{"early_*"} integration assay when \code{modality} is
#'   \code{NULL} (joint models).
#' @param assay Assay index or name within \code{experiment}. Default is 1.
#' @param colors Named character vector mapping risk groups to colours. Default
#'   is \code{c("High risk" = "red", "Low risk" = "deepskyblue")}.
#' @param rho Numeric. For \code{model = "AFTCoop"} only: selects the
#'   \code{beta_rho_<rho>} feature set. If \code{NULL} (default) and only one
#'   rho was trained, it is resolved automatically.
#' @param scale Logical. If \code{TRUE}, features are standardised and
#'   zero-variance features are dropped. Default is \code{TRUE}.
#' @param n_neighbors Integer UMAP neighbourhood size. Default is 15.
#' @param seed Integer random seed for reproducibility. Default is 42.
#' @param ... Additional arguments passed to \code{\link[uwot]{umap}}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{patient_stratification}}, \code{\link{plot_km}},
#'   \code{\link{top_features}}, \code{\link{selected_features}}
#'
#' @export
plot_umap <- function(mae, x, model,
                      modality = NULL,
                      top_feat = NULL,
                      mae_train = NULL,
                      features = NULL,
                      experiment = NULL,
                      assay = 1,
                      colors = .risk_colors(),
                      rho = NULL,
                      scale = TRUE,
                      n_neighbors = 15,
                      seed = 42,
                      ...) {
  strat <- .as_strata(x)
  .check_strata(strat)
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for plot_umap().", call. = FALSE)
  }
  # Project on the assay matching the modality; joint models (modality NULL)
  # default to the integrated assay where their prefixed features live.
  if (is.null(experiment)) experiment <- .default_experiment(mae, modality)
  src <- if (is.null(mae_train)) mae else mae_train

  features <- .resolve_model_features(src, model, modality, top_feat, features,
                                      rho = rho)
  fm       <- .feature_matrix(mae, strat, features, experiment, assay)
  X        <- t(as.matrix(fm$mat))                 # samples x features
  grp      <- fm$groups
  if (nrow(X) < 3) {
    stop("Too few aligned samples (", nrow(X), ") for UMAP.", call. = FALSE)
  }

  X[is.na(X)] <- 0                                 # minimal imputation
  if (scale) {
    v <- apply(X, 2, stats::var)
    X <- X[, is.finite(v) & v > 0, drop = FALSE]
    X <- scale(X)
  }

  set.seed(seed)
  emb <- uwot::umap(X, n_neighbors = min(n_neighbors, nrow(X) - 1L), ...)
  d   <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2], risk_group = grp)

  ggplot2::ggplot(
    d, ggplot2::aes(.data$UMAP1, .data$UMAP2, colour = .data$risk_group)) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = colors, name = NULL) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_bw()
}

# Models whose stored prediction is an AFT-style linear predictor: higher =
# longer survival = LOWER risk. Their score must be negated before thresholding
# so that "higher score = higher risk" holds, matching the Cox models.
#' @keywords internal
#' @noRd
.is_aft_model <- function(model) {
  model %in% c("aft_pen", "AFTCoop")
}

# For AFTCoop, pred and selected_features are lists keyed "beta_rho_<value>".
# Extract the single vector for the requested rho (or error clearly).
#' @keywords internal
#' @noRd
.resolve_aftcoop_rho <- function(x, rho, what = "pred") {
  if (!is.list(x)) return(x)                        # not multi-rho: pass through
  if (is.null(rho)) {
    # Auto-select the single rho present; ask if ambiguous.
    if (length(x) == 1L) return(x[[1L]])
    stop("AFTCoop stores one ", what, " per rho (", paste(names(x), collapse=", "),
         "). Specify `rho` (e.g. rho = 0.5).", call. = FALSE)
  }
  key <- paste0("beta_rho_", rho)
  if (!key %in% names(x))
    stop("rho key '", key, "' not found in AFTCoop ", what, ". Available: ",
         paste(names(x), collapse = ", "), call. = FALSE)
  x[[key]]
}

# benchmark_models() labels AFTCoop rows as "AFTCoop_beta_rho_<v>" when
# aftcoop_rho is NULL. Map those back to "AFTCoop" so downstream readers
# (which look up "AFTCoop" in metadata) can find the model.
#' @keywords internal
#' @noRd
.normalise_model_name <- function(model) {
  if (startsWith(model, "AFTCoop_")) "AFTCoop" else model
}


#' @keywords internal
#' @noRd
.load_orgdb <- function(org_db) {
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Package '", org_db, "' is required for enrichment.", call. = FALSE)
  }
  getExportedValue(org_db, org_db)
}

# Map gene symbols (or another key type) to ENTREZ IDs via clusterProfiler.
#' @keywords internal
#' @noRd
.to_entrez <- function(genes, key_type, orgdb) {
  if (identical(key_type, "ENTREZID")) return(unique(stats::na.omit(genes)))
  map <- suppressWarnings(
    clusterProfiler::bitr(genes, fromType = key_type,
                          toType = "ENTREZID", OrgDb = orgdb))
  unique(map$ENTREZID)
}

# ---------------------------------------------------------------------------
# 1. enrich_features
# ---------------------------------------------------------------------------

#' Functional Enrichment of a Model's Features
#'
#' Runs over-representation analysis (ORA) on the genes of a survival model:
#' Gene Ontology and KEGG via \pkg{clusterProfiler}, Reactome via
#' \pkg{ReactomePA}. Features are obtained exactly as in \code{\link{plot_umap}}
#' / \code{\link{plot_heatmap}} — \code{selected_features()} for penalised
#' models and \code{top_features()} for importance-based models
#' (\code{"blockForest"}, \code{"flexynesis"}) — or supplied directly via
#' \code{features}. Gene symbols are mapped to ENTREZ IDs for KEGG and Reactome;
#' GO is queried on the original key type.
#'
#' @param mae A \code{MultiAssayExperiment} whose \code{metadata()$analysis}
#'   holds the fitted model (typically the training object).
#' @param model Character string naming the model (e.g. \code{"cox_elastic_net"},
#'   \code{"blockForest"}, \code{"flexynesis"}).
#' @param modality Character string naming the omics modality used to locate a
#'   penalised model's selected features (ignored for \code{"blockForest"} /
#'   \code{"flexynesis"}). Default is \code{NULL}: the single per-run
#'   integration tag (\code{"early_<omics>"}) is resolved automatically; pass
#'   the tag (or an omic name) explicitly to disambiguate.
#' @param top_feat Integer number of top features. Used \emph{only} for
#'   \code{"blockForest"} / \code{"flexynesis"}. If \code{NULL} (default),
#'   \code{\link{top_features}}'s own default is used.
#' @param features Optional character vector of gene identifiers, overriding the
#'   automatic retrieval from \code{model}.
#' @param rho Numeric. For \code{model = "AFTCoop"} only: selects the
#'   \code{beta_rho_<rho>} feature set. If \code{NULL} (default) and only one
#'   rho was trained, it is resolved automatically.
#' @param strip_prefix Optional regular expression removed from the start of
#'   each feature name before mapping (e.g. \code{"^(mRNA|meth)_"} when assay
#'   rows are layer-prefixed). Default is \code{NULL}.
#' @param databases Which analyses to run; any of \code{"GO"}, \code{"KEGG"},
#'   \code{"Reactome"}. Default is all three.
#' @param ont GO ontology: \code{"BP"}, \code{"MF"}, \code{"CC"} or
#'   \code{"ALL"}. Default is \code{"ALL"}.
#' @param key_type Key type of the input features for GO and ID mapping. Default
#'   is \code{"SYMBOL"}.
#' @param org_db Name of the organism annotation package. Default is
#'   \code{"org.Hs.eg.db"}.
#' @param kegg_organism KEGG organism code. Default is \code{"hsa"}.
#' @param reactome_organism ReactomePA organism. Default is \code{"human"}.
#' @param pvalue_cutoff,qvalue_cutoff Significance cut-offs passed to each
#'   enrichment call. Defaults are 0.05 and 0.2.
#' @param ... Additional arguments forwarded to the enrichment functions.
#'
#' @return A named \code{list} of \code{enrichResult} objects (one per requested
#'   database), plus \code{genes} (the cleaned input), \code{entrez} (the mapped
#'   IDs) and \code{universe} (the background used, or \code{NULL}). KEGG and
#'   Reactome results are made \code{readable} (ENTREZ back to symbols).
#'
#' @seealso \code{\link{plot_enrichment}}, \code{\link{selected_features}},
#'   \code{\link{top_features}}
#'
#' @export
enrich_features <- function(mae, model,
                            modality = NULL,
                            top_feat = NULL,
                            features = NULL,
                            rho = NULL,
                            strip_prefix = NULL,
                            databases = c("GO", "KEGG", "Reactome"),
                            ont = "ALL",
                            key_type = "SYMBOL",
                            org_db = "org.Hs.eg.db",
                            kegg_organism = "hsa",
                            reactome_organism = "human",
                            pvalue_cutoff = 0.05,
                            qvalue_cutoff = 0.2,
                            ...) {
  databases <- match.arg(databases, several.ok = TRUE)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required for enrich_features().",
         call. = FALSE)
  }
  orgdb <- .load_orgdb(org_db)

  genes <- .resolve_model_features(mae, model, modality, top_feat, features,
                                   rho = rho)
  if (!is.null(strip_prefix)) genes <- sub(strip_prefix, "", genes)
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (!length(genes)) stop("No features to enrich.", call. = FALSE)

  need_entrez <- any(c("KEGG", "Reactome") %in% databases)
  entrez <- if (need_entrez) .to_entrez(genes, key_type, orgdb) else NULL

  out <- list()

  if ("GO" %in% databases) {
    out$GO <- clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = orgdb, keyType = key_type, ont = ont,
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
      readable = TRUE, ...)
  }

  if ("KEGG" %in% databases) {
    kk <- clusterProfiler::enrichKEGG(
      gene = entrez, organism = kegg_organism,
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff, ...)
    if (!is.null(kk)) {
      kk <- clusterProfiler::setReadable(kk, orgdb, keyType = "ENTREZID")
    }
    out$KEGG <- kk
  }

  if ("Reactome" %in% databases) {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      stop("Package 'ReactomePA' is required for Reactome enrichment.",
           call. = FALSE)
    }
    rr <- ReactomePA::enrichPathway(
      gene = entrez, organism = reactome_organism,
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
      readable = TRUE, ...)
    out$Reactome <- rr
  }

  out$genes    <- genes
  out$entrez   <- entrez
  out
}

# ---------------------------------------------------------------------------
# 2. plot_enrichment
# ---------------------------------------------------------------------------

#' Dotplot of an Enrichment Result
#'
#' Convenience wrapper around \code{\link[enrichplot]{dotplot}} for one database
#' of an \code{\link{enrich_features}} result.
#'
#' @param enrichment The list returned by \code{\link{enrich_features}}, or a
#'   single \code{enrichResult} object.
#' @param database Which result to plot when \code{enrichment} is a list:
#'   \code{"GO"}, \code{"KEGG"} or \code{"Reactome"}. Default is \code{"GO"}.
#' @param show_category Integer number of top terms to display. Default is 15.
#' @param ... Additional arguments passed to \code{\link[enrichplot]{dotplot}}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{enrich_features}}
#'
#' @export
plot_enrichment <- function(enrichment,
                            database = c("GO", "KEGG", "Reactome"),
                            show_category = 15,
                            ...) {
  database <- match.arg(database)
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required for plot_enrichment().",
         call. = FALSE)
  }

  res <- if (methods::is(enrichment, "enrichResult")) enrichment
  else enrichment[[database]]

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    stop("No enriched terms to plot for '", database, "'.", call. = FALSE)
  }
  enrichplot::dotplot(res, showCategory = show_category, ...) +
    ggplot2::ggtitle(database)
}

#' @importFrom stats IQR coef mad optim p.adjust predict quantile runif setNames
#' @importFrom utils head type.convert
#' @importFrom graphics axis
#' @importFrom survival coxph basehaz concordance
NULL

# ---------------------------------------------------------------------------
# 3. plot_heatmap
# ---------------------------------------------------------------------------

#' Feature Heatmap Split by Risk Group
#'
#' Draws a \code{\link[pheatmap]{pheatmap}} of a model's features across the
#' stratified patients. Features are selected exactly as in
#' \code{\link{plot_umap}}: \code{selected_features()} for penalised models and
#' \code{top_features()} for importance-based models
#' (\code{"blockForest"}, \code{"flexynesis"}), chosen from \code{model}.
#'
#' By default the columns (patients) are clustered, so the risk-group split is
#' \emph{not} imposed: the \code{Risk} annotation bar then shows whether the
#' expression-based clusters line up with the risk groups. Set
#' \code{cluster_cols = FALSE} to instead order patients by risk group
#' (\code{"High risk"} first) with a column gap marking the split. Rows
#' (features) are clustered with Ward's method. Values are row z-scored and
#' winsorised to \code{[-cap, cap]} so a few outliers do not wash out the
#' palette. Feature information is read from \code{mae_train} (a training
#' artifact); the matrix shown is drawn from \code{mae}.
#'
#' @param mae A \code{MultiAssayExperiment} providing the expression matrix
#'   (typically the stratified test object).
#' @param x Either the \code{list} returned by
#'   \code{\link{patient_stratification}} or its \code{strata}
#'   \code{data.frame}.
#' @param model Character string naming the model whose features are shown.
#' @param modality Character string naming the omics modality. Locates a
#'   penalised model's selected features (ignored for \code{"blockForest"} /
#'   \code{"flexynesis"}) and, unless \code{experiment} is given, selects the
#'   assay to draw. Default is \code{NULL}: the single per-run integration tag
#'   (\code{"early_<omics>"}) is resolved automatically; pass the tag (or an
#'   omic name) explicitly to disambiguate.
#' @param top_feat Integer number of top features. Used \emph{only} for
#'   \code{"blockForest"} / \code{"flexynesis"}. If \code{NULL} (default),
#'   \code{\link{top_features}}'s own default is used.
#' @param mae_train A \code{MultiAssayExperiment} from which to read the model's
#'   features. Defaults to \code{mae}.
#' @param features Optional character vector of feature names, overriding the
#'   automatic retrieval from \code{model}.
#' @param rho Numeric. For \code{model = "AFTCoop"} only: selects the
#'   \code{beta_rho_<rho>} feature set. If \code{NULL} (default) and only one
#'   rho was trained, it is resolved automatically.
#' @param experiment Name or index of the experiment to draw. If \code{NULL}
#'   (default), the assay matching \code{modality} is used, or the single
#'   \code{"early_*"} integration assay when \code{modality} is \code{NULL}
#'   (joint models).
#' @param assay Assay index or name within \code{experiment}. Default is 1.
#' @param colors Named character vector mapping risk groups to annotation
#'   colours. Default is
#'   \code{c("High risk" = "red", "Low risk" = "deepskyblue")}.
#' @param annotations Optional character vector of \code{colData(mae)} column
#'   names to add as extra column annotations alongside \code{Risk} (e.g.
#'   \code{c("status", "pathologic_stage")}). Colours for these are assigned by
#'   \code{pheatmap}.
#' @param discrete_max Integer. Extra annotation columns are shown as discrete
#'   categories unless they are numeric with more than \code{discrete_max}
#'   distinct values, which stay a continuous gradient. Default is 5.
#' @param heat_colors Character vector of colours for the heatmap body. Default
#'   is a navy-white-firebrick ramp of 200 steps.
#' @param cap Numeric. Row z-scores are clipped to \code{[-cap, cap]}. Default
#'   is 3.
#' @param cluster_cols Logical. If \code{TRUE} (default), patients are clustered
#'   and their order is left to \code{pheatmap}; if \code{FALSE}, they are
#'   ordered by risk group with a gap marking the split.
#' @param cluster_rows Logical, whether to cluster features. Default is
#'   \code{TRUE}.
#' @param clustering_method Agglomeration method for row clustering. Default is
#'   \code{"ward.D2"}.
#' @param show_rownames,show_colnames Logical, whether to print feature / sample
#'   labels. Defaults are \code{TRUE} and \code{FALSE}.
#' @param fontsize_row Font size for feature labels. Default is 6.
#' @param ... Additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return The \code{pheatmap} object (invisibly), drawn as a side effect.
#'
#' @seealso \code{\link{patient_stratification}}, \code{\link{plot_umap}},
#'   \code{\link{top_features}}, \code{\link{selected_features}}
#'
#' @export
plot_heatmap <- function(mae, x, model,
                         modality = NULL,
                         top_feat = NULL,
                         mae_train = NULL,
                         features = NULL,
                         rho = NULL,
                         experiment = NULL,
                         assay = 1,
                         colors = .risk_colors(),
                         annotations = NULL,
                         discrete_max = 5,
                         heat_colors = grDevices::colorRampPalette(
                           c("navy", "white", "firebrick3"))(200),
                         cap = 3,
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         clustering_method = "ward.D2",
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         fontsize_row = 6,
                         ...) {
  strat <- .as_strata(x)
  .check_strata(strat)
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for plot_heatmap().", call. = FALSE)
  }
  if (is.null(experiment)) experiment <- .default_experiment(mae, modality)
  src <- if (is.null(mae_train)) mae else mae_train

  features <- .resolve_model_features(src, model, modality, top_feat, features,
                                      rho = rho)
  fm       <- .feature_matrix(mae, strat, features, experiment, assay)

  # When clustering columns, let pheatmap order the samples; otherwise order
  # them by risk group (High risk first) and mark the split with a gap.
  if (cluster_cols) {
    m   <- as.matrix(fm$mat)
    grp <- fm$groups
    gap <- NULL
  } else {
    ord <- order(fm$groups, decreasing = TRUE)     # levels: High risk, Low risk
    m   <- as.matrix(fm$mat[, ord, drop = FALSE])
    grp <- fm$groups[ord]
    gap <- sum(grp == "High risk")
  }

  # Row z-score, then winsorise to [-cap, cap] with symmetric colour breaks.
  m <- t(scale(t(m)))
  m[is.na(m)] <- 0                                 # zero-variance rows
  m[m >  cap] <-  cap
  m[m < -cap] <- -cap
  breaks <- seq(-cap, cap, length.out = length(heat_colors) + 1L)

  ann     <- data.frame(Risk = grp, row.names = colnames(m))
  if (!is.null(annotations)) {
    cd   <- as.data.frame(MultiAssayExperiment::colData(mae),
                          stringsAsFactors = FALSE)
    miss <- setdiff(annotations, colnames(cd))
    if (length(miss)) {
      stop("Annotation column(s) not found in colData: ",
           paste(miss, collapse = ", "), call. = FALSE)
    }
    ann <- cbind(ann, cd[colnames(m), annotations, drop = FALSE])
    # Categorical unless numeric with many distinct values (keep those as a
    # continuous gradient); this avoids e.g. status 0/1 becoming a gradient.
    for (a in annotations) {
      col <- ann[[a]]
      if (!is.numeric(col) ||
          length(unique(stats::na.omit(col))) <= discrete_max) {
        ann[[a]] <- factor(col)
      }
    }
  }
  ann_col <- list(Risk = colors[levels(grp)])

  ph <- pheatmap::pheatmap(
    m,
    color             = heat_colors,
    breaks            = breaks,
    annotation_col    = ann,
    annotation_colors = ann_col,
    cluster_cols      = cluster_cols,
    cluster_rows      = cluster_rows,
    clustering_method = clustering_method,
    gaps_col          = gap,
    scale             = "none",
    show_rownames     = show_rownames,
    show_colnames     = show_colnames,
    fontsize_row      = fontsize_row,
    border_color      = NA,
    main              = paste0(model, " \u2014 ", modality),
    ...
  )
  invisible(ph)
}
