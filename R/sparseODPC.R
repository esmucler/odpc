#' @title Automatic Choice of Regularization Parameters for Sparse One-Sided Dynamic Principal Components via Cross-Validation
#' @param Z Data matrix. Each column is a different time series.
#' @param h Forecast horizon.
#' @param k_max Maximum possible number of lags to use.
#' @param max_num_comp Maximum possible number of components to compute.
#' @param window_size The size of the rolling window used to estimate the forecasting error.
#' @param method A string specifying the algorithm used. Options are 'ALS', 'mix' or 'gradient'. See details in \code{\link{odpc}}.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.
#' @param ncores Number of cores to parallelise over.



#' @return 
#' An object of class odpcs, that is, a list of length equal to the number of computed components, each computed using the optimal value of k. 
#' The i-th entry of this list is an object of class \code{odpc}, that is, a list with entries
#'\item{f}{Coordinates of the i-th dynamic principal component corresponding to the periods \eqn{k_1 + 1,\dots,T}.}
#'\item{mse}{Mean squared error of the reconstruction using the first i components.}
#'\item{k1}{Number of lags used to define the i-th dynamic principal component f.}
#'\item{k2}{Number of lags of f used to reconstruct.}
#'\item{alpha}{Vector of intercepts corresponding to f.}
#'\item{a}{Vector that defines the i-th dynamic principal component}
#'\item{B}{Matrix of loadings corresponding to f. Row number \eqn{k} is the vector of \eqn{k-1} lag loadings.}
#'\item{call}{The matched call.}
#'\item{conv}{Logical. Did the iterations converge?}
#'\item{lambda}{Regularization parameter used for this componen}
#'\code{components}, \code{fitted}, \code{plot} and \code{print} methods are available for this class.
#' 
#' 
#' @description
#' Computes Sparse One-Sided Dynamic Principal Components, choosing the number of components and regularization parameters automatically, to minimize an estimate
#' of the forecasting mean squared error.
#'
#' @details 
#' TBA
#' 
#' @references
#' Peña D., Smucler E. and Yohai V.J. (2017). “Forecasting Multiple Time Series with One-Sided Dynamic Principal Components.” Available at https://arxiv.org/abs/1708.04705.
#'
#' @seealso \code{\link{odpc}}, \code{\link{crit.odpc}}, \code{\link{forecast.odpcs}}
#' 
#' @examples 
#' T <- 50 #length of series
#' m <- 10 #number of series
#' set.seed(1234)
#' f <- rnorm(T + 1)
#' x <- matrix(0, T, m)
#' u <- matrix(rnorm(T * m), T, m)
#' for (i in 1:m) {
#'   x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
#' }
#' # Choose parameters to perform a one step ahead forecast. Use 1 or 2 lags, only one component 
#' # and a window size of 2 (artificially small to keep computation time low). Use two cores for the
#' # loop over k, two cores for the loop over the window
#' fit <- cv.sparse_odpc(x, h=1, k_max = 1, max_num_comp = 1, window_size = 2, ncores = 1)
cv.sparse_odpc <- function(Z, h, k_max = 3, max_num_comp = 2, window_size, method, tol = 1e-04, niter_max = 500, ncores=1) {
  
  if (all(!inherits(Z, "matrix"), !inherits(Z, "mts"),
          !inherits(Z, "xts"), !inherits(Z, "data.frame"),
          !inherits(Z, "zoo"))) {
    stop("Z should belong to one of the following classes: matrix, data.frame, mts, xts, zoo")
  } else if (any(dim(Z)[2] < 2, dim(Z)[1] < 10)) {
    stop("Z should have at least ten rows and two columns")
  } else if (any(anyNA(Z), any(is.nan(Z)), any(is.infinite(Z)))) {
    stop("Z should not have missing, infinite or nan values")
  }
  if (!inherits(tol, "numeric")) {
    stop("tol should be numeric")
  } else if (!all(tol < 1, tol > 0)) {
    stop("tol be between 0 and 1")
  }
  if (!inherits(niter_max, "numeric")) {
    stop("niter_max should be numeric")
  } else if (any(!niter_max == floor(niter_max), niter_max <= 0)) {
    stop("niter_max should be a positive integer")
  }
  if (!missing(method)) {
    if (!(method %in% c('ALS', 'mix', 'gradient'))) {
      stop("method should be ALS, mix or gradient")
    }
  }
  
  if (missing(method)){
    if (ncol(Z) > 10){
      method_num <- 3 # Use gradient method for moderately fat data sets
    } else {
      method_num <- 1
    }
  } else {
    method_num <- switch(method, 'ALS' = 1, 'mix' = 2, 'gradient' = 3)  
  }
  if (missing(window_size)){
    window_size <- floor(0.2 * nrow(Z))
  }
  cl <- makeCluster(getOption("cl.cores", ncores))
  old_best_mse <- Inf # Previous estimate of best forecast MSE
  
  Z_train <- Z[1:(nrow(Z) - window_size),]
  response <- Z_train[(2 * k_max + 1):(nrow(Z_train)),]
  
  odpc_fit <- odpc(Z=Z_train, ks=k_max, method=method, tol=tol, niter_max=niter_max)
  sparse_path <- sparse_odpc_path(fit_component=odpc_fit[[1]], Z=Z_train, response=response)
  sparse_path <- lapply(sparse_path, function(fit) {list(fit)})
  
  forecasts <- forecast_sparse_rolled(sparse_path=sparse_path, Z=Z, window_size=window_size, h=h, cl=cl)
  best_fit <- get_best_sparse_fit(forecasts=forecasts, sparse_path=sparse_path, Z=Z, h=h)
  
  opt_comp <- best_fit$opt_comp
  new_best_mse <- best_fit$opt_mse
  new_opt_lambda <- opt_comp[[1]]$lambda
  lambdas <- c(new_opt_lambda)
  opt_comp <- construct.odpcs(out=opt_comp, data=Z_train, fn_call=match.call())
  num_comp <- 1
  
  while (new_best_mse < old_best_mse & num_comp < max_num_comp){
    
    num_comp <- num_comp + 1
    response_residual <- response - fitted(opt_comp, num_comp=length(opt_comp))
    
    odpc_fit <- convert_rename_comp(odpc_priv(Z = Z_train, resp=response_residual, k_tot_max=2*k_max,
                                    k1 = k_max, k2 = k_max, num_comp=num_comp, f_ini = 0,
                                    passf_ini = FALSE, tol = tol, niter_max = niter_max, method=method_num), wrap=TRUE)
    odpc_fit <- construct.odpcs(odpc_fit, data=Z_train, fn_call=match.call())[[1]]
    
    sparse_path <- sparse_odpc_path(fit_component=odpc_fit, Z=Z_train, response=response_residual)
    sparse_path <- lapply(sparse_path, function(fit) {out <- append(opt_comp, list(fit)); class(out) <- append(class(out), 'odpcs'); return(out)})
    
    forecasts <- forecast_sparse_rolled(sparse_path=sparse_path, Z=Z, window_size=window_size, h=h, cl=cl)
    best_fit <- get_best_sparse_fit(forecasts=forecasts, sparse_path=sparse_path, Z=Z, h=h)
    
    
    if (best_fit$opt_mse < new_best_mse){
      opt_comp <- best_fit$opt_comp
      new_best_mse <- best_fit$opt_mse
      new_opt_lambda <- opt_comp[[num_comp]]$lambda
      lambdas <- c(lambdas, new_opt_lambda)
      opt_comp <- construct.odpcs(out=opt_comp, data=Z_train, fn_call=match.call())
    } else {
      old_best_mse <- new_best_mse
    }
  }
  
  response_full <- Z[(2*k_max + 1):nrow(Z),]
  final_sparse_fit <- vector(length=length(lambdas), mode='list')
  # TODO test the output of this loop by checking that the mses of the reconsutrctions match
  for (iter in 1:length(lambdas)){
    odpc_fit <- convert_rename_comp(odpc_priv(Z = Z, resp=response_full, k_tot_max=2*k_max,
                                         k1 = k_max, k2 = k_max, num_comp=iter, f_ini = 0,
                                         passf_ini = FALSE, tol = tol, niter_max = niter_max, method=method_num), wrap=TRUE)
    odpc_fit <- construct.odpcs(odpc_fit, data=Z, fn_call=match.call())[[1]]
    
    final_sparse_fit[[iter]] <- sparse_odpc_path(fit_component=odpc_fit, Z=Z, response=response_full, lambda=lambdas[iter])[[1]]
    response_full <- response_full - fitted(final_sparse_fit[[iter]])
  }
  final_sparse_fit <- construct.odpcs(out=final_sparse_fit, data=Z, fn_call=match.call())
  on.exit(stopCluster(cl))
  return(final_sparse_fit)
}

get_best_sparse_fit <- function(sparse_path, forecasts, Z, h){
  mses <- get_ave_mse_sparse(forecasts=forecasts, Z=Z, h=h)
  opt_lambda_ind <- which.min(mses)
  opt_comp <- sparse_path[[opt_lambda_ind]]
  return(list(opt_comp=opt_comp, opt_mse=min(mses)))
}

get_ave_mse_sparse <- function(forecasts, Z, h){
  centered_forecasts <- lapply(seq_along(forecasts), function(w) {scale(t(forecasts[[w]]), center=Z[nrow(Z) - w + 1,], scale=FALSE)})
  mses_per_window <- sapply(seq_along(centered_forecasts), function(w) { apply(centered_forecasts[[w]] , 1, function(row) mean(row**2))})
  mses <- apply(mses_per_window, 1, mean)
}

forecast_sparse_rolled <- function(sparse_path, Z, window_size, h, cl){
  N <- nrow(Z)
  forecasts <- parLapply(cl=cl, X=1:window_size, fun=function(w) {forecast_sparse_path(sparse_path=sparse_path, rolled_data=Z[1:(N - h - w + 1),], h=h)})
  return(forecasts)
}

forecast_sparse_path <- function(sparse_path, rolled_data, h){
  # each column is a different lambda
  forecasts <- sapply(sparse_path, function(fit) {forecast_sparse_odpc(fit=fit, rolled_data=rolled_data, h=h)})
  return(forecasts)
}

forecast_sparse_odpc <- function(fit, rolled_data, h){
  ncomp <- length(fit)
  fore <- 0
  # sparse fits have the same k1 and k2 in all components
  k_max <- fit[[1]]$k1
  # get components defined using a in fit but data in rolled_data
  new_comps <- sapply(fit, function(fit_component) {get_new_comp(a=fit_component$a, rolled_data=rolled_data, k_max=k_max)})
  # TODO pass ... args to auto.arima
  fores_comps <- apply(new_comps, 2, function(x, h, ...) { auto <- auto.arima(x)
                                                      return(forecast(auto, h)$mean[h])
                                                      }, h)
  new_comps <- rbind(new_comps, fores_comps)
  for (i in 1:ncomp){
    matF <- getMatrixFore(f=new_comps[, i], k2=k_max, h=h)
    fore <- fore + matF %*% rbind(as.vector(fit[[i]]$alpha), as.matrix(fit[[i]]$B))
  }
  return(fore)
}

get_new_comp <- function(a, rolled_data, k_max){
  # Computes new f using the a passed as input and the data in rolled_data
  matF <- getMatrixF_sparse_forecast(Z=rolled_data, k1=k_max, k2=k_max, k_tot=2*k_max, a=a)
  new_comp <- rep(NA, nrow(rolled_data) - k_max)
  new_comp[1:k_max] <- matF[1:k_max, k_max+2]
  new_comp[(k_max + 1):(nrow(rolled_data) - k_max)] <- matF[,2]
  return(new_comp)
}

sparse_odpc_path <- function(fit_component, Z, response, lambda=NULL, ...){
  N <- nrow(Z)
  k1 <- fit_component$k1
  k2 <- fit_component$k2
  matreg <- getMatrixZj0(Z=Z, k1=k1, k_tot=k1, j=k1)
  reg_path <- glmnet(y=fit_component$f, x=matreg, lambda=lambda, intercept=FALSE, ...)
  lambda_fitted <- reg_path$lambda
  component_path <- predict.glmnet(reg_path, newx=matreg)
  coordinates_path <- coef.glmnet(reg_path)
  matD_path <- apply(component_path, 2, function(comp) {update_loadings(comp, response, k1=k1, k2=k2)})
  sparse_odpc_path <- lapply(seq_along(lambda_fitted), function(k) {build_sparse_component(mse=matD_path[[k]]$mse, matD=matD_path[[k]]$matD, lambda=lambda_fitted[k],
                                                                                component=component_path[,k], a=coordinates_path[,k], fit=fit_component)}) 
  fit_component$lambda <- 0
  sparse_odpc_path <- append(sparse_odpc_path, list(fit_component))
  return(sparse_odpc_path)
}

update_loadings <- function(component, response, k1, k2){
  matF <- getMatrixFitted(f=component, k1=k1, k2=k2)
  matD <- MASS::ginv(matF) %*% response
  mse <- mean((response - matF %*% matD)**2)
  return(list(matD=matD, mse=mse))
}

build_sparse_component <- function(mse, matD, component, a, lambda, fit){
  fit$mse <- mse
  fit$B <- matD[2:nrow(matD),]
  fit$alpha <- matD[1,]
  fit$f <- component
  fit$a <- a[-1]
  fit$lambda <- lambda
  return(fit)
}