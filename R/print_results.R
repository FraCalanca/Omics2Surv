print_single_result = function(train_obj, test_obj, method_list, lambda_type= NULL) {
  result_table = tibble::tibble(

    p = if (method=="penAFT"){
      length(train_obj$coefficients)
    } else{
      nrow(train_obj$final_model$beta)
    },

    Method = method_list,

    Lambda = if (method %in% c("lasso", "EN", "adaptiveLASSO", "penAFT", "CoopCox")) {
      lambda_type
    },

    Alpha = if (method %in% c("EN", "penAFT")) {
      train_obj$best_alpha
    } else {
      1
    },

    Best_Lambda = if (method %in% c("lasso", "EN", "adaptiveLASSO", "penAFT", "CoopCox")) {
      round(train_obj$selected_lambda,5)
    },

    d = if (method %in% c("lasso", "EN", "adaptiveLASSO", "penAFT", "CoopCox", "AFTCoop")) {
      length(train_obj$selected_features)
    },

    CIndex_Training = round(train_obj$train_cindex,4),

    CIndex_Test = round(test_obj$Cindex_test,4)

  )

  print(result_table)

}

print_results_comparison = function(train_list, test_list) {

  if (length(train_list) != length(test_list)) {

    stop("Error: different length for train_list and test_list !")

  }

  methods = c(rep("Lasso", 4), rep("Elastic Net", 4), rep("adaptiveLASSO",4), "penAFT")

  lambda_types = c(rep(c("min", "1se"), 6), "min")

  type_measures = c(rep(c("C", "C", "deviance", "deviance"), 3), "NA")

  results_table = tibble::tibble(


    Method = methods,

    Lambda = lambda_types,

    Measure = type_measures,

    Alpha_Sel = c(

      rep(1, 4),

      sapply(train_list[5:8], function(x) x$best_alpha),

      rep(1, 4),

      sapply(train_list[13], function(x) x$best_alpha)

    ),

    Best_Lambda = sapply(train_list, function(x) round(x$selected_lambda,5)),

    p = c(
      sapply(train_list[1:12], function(x) nrow(x$final_model$beta)),
      sapply(train_list[13], function(x) length(x$coefficients))
    ),

    d = sapply(train_list, function(x) length(x$selected_features)),

    CIndex_Training = sapply(train_list, function(x) round(x$train_cindex,4)),

    CIndex_Test = sapply(test_list, function(x) round(x$Cindex_test,4))

  )

  print(results_table)

}
