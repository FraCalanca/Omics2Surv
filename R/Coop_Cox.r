#' Train Cooperative Learning Model for Multi-Omics Survival Analysis
#'
#' Trains a cooperative Lasso Cox regression model integrating 2 or 3 omics
#' layers. An internal validation subset is held out for lambda selection via
#' cross-validation. Can be called directly or via \code{\link{train_j_models}}.
#'
#' @param X_list Named list of 2 or 3 omics matrices (samples x features).
#'   Row names must be sample identifiers, consistent across all layers.
#' @param y_training Data frame with \code{overall_survival} and \code{status}
#'   columns. Row names must match sample identifiers in \code{X_list}.
#' @param val_fraction Numeric proportion of training samples to hold out as
#'   internal validation for lambda selection. Default is \code{0.1}.
#' @param lambdaVector Numeric vector of lambda values to evaluate.
#'   Default is \code{seq(0.1, 0.9, 0.1)}.
#' @param folds Integer number of cross-validation folds. Default is 5.
#' @param alpha Numeric cooperation parameter: 1 for 2 layers, 0.5 for 3.
#'   Default is 1.
#' @param maxit Integer maximum iterations. Default is 1000.
#' @param seed Integer seed. Default is 123.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{model}: Trained cooplasso model.
#'     \item \code{selected_lambda}: Optimal lambda selected by CV.
#'     \item \code{selected_features}: Non-zero coefficient feature names.
#'     \item \code{pred}: Named numeric vector of risk scores for the full
#'       input set (train + internal validation reassembled in original order).
#'     \item \code{omics}: Names of omics layers used.
#'   }
#'
#' @seealso \code{\link{test_cooplearning}}, \code{\link{train_j_models}}
#'
#' @references
#' Hahn G, Prokopenko D, Hecker J, Lutz SM, Mullin K, Sejour L, Hide W,
#' Vlachos I, DeSantis S, Tanzi RE, Lange C. Prediction of disease-free
#' survival for precision medicine using cooperative learning on multi-omic
#' data. \emph{Briefings in Bioinformatics}, 25(4):bbae267, July 2024.
#' \doi{10.1093/bib/bbae267}
#'
#' @export
train_cooplearning <- function(
    X_list,
    y_training,
    val_fraction = 0.1,
    lambdaVector = seq(0.1, 0.9, 0.1),
    folds        = 5,
    alpha        = 1,
    maxit        = 100,
    seed         = 123
) {

  set.seed(seed)

  stopifnot(is.list(X_list))
  if (!length(X_list) %in% c(2L, 3L)) {
    stop("'X_list' must contain exactly 2 or 3 omics layers.")
  }

  # ── Align samples across layers and y ────────────────────────────
  common_ids <- Reduce(intersect,
                       c(lapply(X_list, rownames), list(rownames(y_training))))
  if (length(common_ids) == 0L) {
    stop("No common samples across all omics layers and y_training.")
  }
  X_list     <- lapply(X_list, function(x) x[common_ids, , drop = FALSE])
  y_training <- y_training[common_ids, , drop = FALSE]

  # ── Internal train / validation split ────────────────────────────
  n         <- length(common_ids)
  val_idx   <- sample(seq_len(n), size = round(n * val_fraction))
  train_idx <- setdiff(seq_len(n), val_idx)

  train_ids <- common_ids[train_idx]
  val_ids   <- common_ids[val_idx]

  X_train_list <- lapply(X_list, function(x) x[train_idx, , drop = FALSE])
  X_val_list   <- lapply(X_list, function(x) x[val_idx,   , drop = FALSE])
  y_train_int  <- y_training[train_idx, , drop = FALSE]

  time  <- as.numeric(y_train_int$overall_survival)
  event <- as.numeric(y_train_int$status)

  Z1  <- X_train_list[[1]]
  Z2  <- X_train_list[[2]]
  Z3  <- if (length(X_train_list) == 3L) X_train_list[[3]] else NULL
  Z_val <- do.call(cbind, X_val_list)

  # ── Cross-validation for lambda ───────────────────────────────────
  message("Cross-validation for CoopCox optimal lambda...")
  opt_lambda <- crossval(
    T            = time,
    y            = event,
    X1           = Z1,
    X2           = Z2,
    X3           = Z3,
    lambdaVector = lambdaVector,
    folds        = folds,
    alpha        = alpha,
    maxit        = maxit
  )
  message("Optimal lambda: ", opt_lambda)

  # ── Final model ───────────────────────────────────────────────────
  message("Training CoopCox...")
  final_model <- cooplasso(
    Z1         = Z1,
    Z2         = Z2,
    Z3         = Z3,
    T          = time,
    delta      = event,
    validation = Z_val,
    lambda     = opt_lambda,
    alpha      = alpha,
    maxit      = maxit,
    normalize  = TRUE
  )

  # ── Predictions reassembled in original order ─────────────────────
  X_train_combined <- do.call(cbind, X_train_list)
  X_val_combined   <- do.call(cbind, X_val_list)

  pred_train <- setNames(as.numeric(X_train_combined %*% final_model$b), train_ids)
  pred_val   <- setNames(as.numeric(X_val_combined   %*% final_model$b), val_ids)
  pred_full  <- c(pred_train, pred_val)[common_ids]

  # ── Selected features ─────────────────────────────────────────────
  coef_vec          <- setNames(as.numeric(final_model$b), colnames(X_train_combined))
  selected_features <- names(coef_vec)[coef_vec != 0]

  list(
    model             = final_model,
    selected_lambda   = opt_lambda,
    selected_features = selected_features,
    pred              = pred_full,
    omics             = names(X_list)
  )
}


#' Test Cooperative Learning Model
#'
#' Applies a trained CoopCox model to new data. Can be called directly or
#' via \code{\link{predict_j_models}}.
#'
#' Omics layers and features are assumed to be in the same order as in
#' training.
#'
#' @param X_list_test Named list of omics matrices (samples x features) for
#'   the test set. Must have the same omics layers and features, in the same
#'   order, as training.
#' @param training_result List returned by \code{\link{train_cooplearning}}.
#'
#' @return A list with \code{pred}: named numeric vector of risk scores.
#'
#' @seealso \code{\link{train_cooplearning}}, \code{\link{predict_j_models}}
#'
#' @export
test_cooplearning <- function(X_list_test, training_result) {

  stopifnot(is.list(X_list_test))

  X_test_combined <- do.call(cbind, X_list_test)

  pred_test <- setNames(
    as.numeric(X_test_combined %*% training_result$model$b),
    rownames(X_test_combined)
  )

  list(pred = pred_test)
}

require(glmnet)
require(ROCR)
require(survival)
require(data.table)


# *******************************************
# likelihood functions used by cooplasso/f2
# *******************************************

# vector b (of length ncol(Z)), vector of times T, vector of censoring indicators delta, matrix Z with one row per subject
l1 <- function(b,T,delta,Z) {
  # indicator T matrix, where each row i contains the indicators T>=T[i]
  indicator <- matrix(T,byrow=TRUE,nrow=length(T),ncol=length(T))>=T
  # multiply each row with exp(b*Z[j,]) for j=1,2,...
  temp <- t( t(indicator) * exp(as.numeric(Z %*% b)) )
  res <- delta * (as.numeric(Z %*% b) - log(rowSums(temp)))
  return(sum(res))
}

# numerical integration
integral <- function(x,y) {
  res <- numeric(length(x))
  for(i in 2:length(x)) {
    res[i] <- (x[i]-x[i-1])*min(y[i],y[i-1]) + 0.5*(x[i]-x[i-1])*abs(y[i]-y[i-1])
  }
  return(res)
}

# normalize vector
vectornormalize <- function(x,option=1) {
  if(option==1) return( x/sqrt(sum(x**2)) )
  if(option==2) return( (x-min(x))/(max(x)-min(x)) )
}

# numerical differentiation to compute basehaz from integrated basehaz
numdiff <- function(m) {
  y <- m[,1]
  x <- m[,2]
  y <- diff(y)#/diff(x)
  x <- x[-length(x)]
  cbind(y,x)
}


# ***************************************************
# cooperative learning
# ***************************************************

# likelihood of the cooperative cox-lasso
f2 <- function(b,T,delta,Z1,Z2,Z3,lambda,b_hat,alpha) {
  n <- length(T)
  indices1 <- 1:ncol(Z1)
  b1 <- b[indices1]
  indices2 <- (ncol(Z1)+1):(ncol(Z1)+ncol(Z2))
  b2 <- b[indices2]
  if(is.null(Z3)) {
    return( -1/n*l1(b1,T,delta,Z1) -1/n*l1(b2,T,delta,Z2) + lambda*sum(abs(b)/abs(b_hat)) + alpha*sum((Z1 %*% b1 - Z2 %*% b2)**2) )
  }
  else {
    indices3 <- (ncol(Z1)+ncol(Z2)+1):(ncol(Z1)+ncol(Z2)+ncol(Z3))
    b3 <- b[indices3]
    return( -1/n*l1(b1,T,delta,Z1) -1/n*l1(b2,T,delta,Z2) -1/n*l1(b3,T,delta,Z3) + lambda*sum(abs(b)/abs(b_hat)) + alpha*sum((Z1 %*% b1 - Z2 %*% b2)**2) + alpha*sum((Z2 %*% b2 - Z3 %*% b3)**2) )
  }
}

# full cooperative cox-lasso
# Reference:
#   Hahn G, Prokopenko D, Hecker J, Lutz SM, Mullin K, Sejour L, Hide W,
#   Vlachos I, DeSantis S, Tanzi RE, Lange C. Prediction of disease-free
#   survival for precision medicine using cooperative learning on multi-omic
#   data. Briefings in Bioinformatics, 25(4):bbae267, July 2024.
#   https://doi.org/10.1093/bib/bbae267
cooplasso <- function(Z1,Z2,Z3=NULL,T,validation,delta,lambda,alpha,maxit=1000,normalize=TRUE) {
  Z <- cbind(Z1,Z2,Z3)
  b_hat <- optim(par=runif(ncol(Z)), fn=function(b) -l1(b,T,delta,Z), method="SANN", gr=NULL, control=list(maxit=maxit/10))$par
  if(normalize) b_hat <- vectornormalize(b_hat)
  b <- optim(par=runif(ncol(Z)), fn=function(b) f2(b,T,delta,Z1,Z2,Z3,lambda,b_hat,alpha), method="SANN", gr=NULL, control=list(maxit=maxit))$par
  if(normalize) b <- vectornormalize(b)
  # compute baseline hazard and multiply
  coxfit <- coxph(Surv(T,delta==1)~0)
  bhest <- numdiff(basehaz(coxfit))
  xv <- bhest[,2]
  hv <- sapply(1:nrow(validation), function(i) bhest[,1]*exp(sum(b*validation[i,])) )
  # calculate survival
  sv <- sapply(1:ncol(hv), function(i) exp(-cumsum(integral(xv,hv[,i]))) )
  return(list(xv=xv,hv=hv,sv=sv, b=b))
}


# *******************************************************
# cross validation to determine suitable lambda for Lasso
# *******************************************************

# find closest index to given value in vector
closestIndex <- function(x,value) {
  which.min(abs(x-value))[1]
}

crossval <- function(T,y,X1,X2,X3=NULL,lambdaVector=seq(0.1,1,0.1),folds=10,maxit=100,alpha=0.5,makeplot=FALSE) {
  resL1 <- matrix(0,nrow=folds,ncol=length(lambdaVector))
  resAUC <- matrix(0,nrow=folds,ncol=length(lambdaVector))
  resC <- matrix(0,nrow=folds,ncol=length(lambdaVector))
  resIBS <- matrix(0,nrow=folds,ncol=length(lambdaVector))
  pool <- sample(1:nrow(X1))
  size <- nrow(X1) %/% folds
  for(f in 1:folds) {
    if(makeplot) print(c(f,folds))
    m <- min(size,length(pool))
    indices <- pool[1:m]
    pool <- pool[-(1:m)]

    Z1 <- X1[indices,]
    Z2 <- X2[indices,]
    Z3 <- X3[indices,]
    Ztrain <- cbind(Z1,Z2,Z3)
    Ttrain <- T[indices]
    Zvalid <- cbind(X1[-indices,],X2[-indices,],X3[-indices,])
    Tvalid <- T[-indices]
    # delta is I(T<=C)
    delta <- y[indices]
    status <- y[-indices]

    # loop
    for(i in 1:length(lambdaVector)) {
      lambda <- lambdaVector[i]
      if(makeplot) print(lambda)
      res <- cooplasso(Z1,Z2,Z3,Ttrain,Zvalid,delta,lambda=lambda,alpha=alpha,maxit=maxit)

      # get survival prediction
      temp <- numeric(nrow(Zvalid))
      for(j in 1:nrow(Zvalid)) {
        index <- closestIndex(res$xv,Tvalid[j])
        temp[j] <- res$sv[index,j]
      }

      # metrics
      indices <- is.na(temp)
      temp[indices] <- 0
      resL1[f,i] <- sum(abs(status-temp))/length(status)
      p <- prediction(temp,status)
      resAUC[f,i] <- performance(p,"auc")@y.values[[1]]
      resC[f,i] <- concordance(y~x,data=data.frame(y=status,x=temp))$concordance
      resIBS[f,i] <- sum((status-temp)**2)/length(status)
    }
  }
  metric <- colMeans(resIBS)
  if(makeplot) {
    plot(lambdaVector,metric,type="p",pch="+",xlab="lambda",ylab="Integrated Brier Score",xaxt="n",xlim=c(0,1))
    axis(side=1,at=(0:10)/10,labels=(0:10)/10)
  }
  else {
    return(lambdaVector[which.min(metric)])
  }
}
