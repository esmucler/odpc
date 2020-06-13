#' @title Automatic Choice of Regularization Parameters for Sparse One-Sided Dynamic Principal Components via Cross-Validation
#' @param Z Data matrix. Each column is a different time series.
#' @param h Forecast horizon.
#' @param k_list List of values of k to choose from.
#' @param nlambda Length of penalty sequence.
#' @param alpha_en Between 0 and 1, elastic net
#' @param window_size The size of the rolling window used to estimate the forecasting error.
#' @param ks Optional, number of lags to use.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.
#' @param eps Between 0 and 1, used to build penalty sequence
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
#' fit <- cv.sparse_odpc(x, h=1, k_list = 1, window_size = 2, ncores = 1)
cv.sparse_odpc <- function(Z, h, k_list = 1:3, nlambda=20, alpha_en=0.95, window_size, ks, tol = 1e-04, niter_max = 500, eps=1e-3, ncores=1) {
  
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
  
  if (missing(window_size)){
    window_size <- floor(0.2 * nrow(Z))
  }
  cl <- makeCluster(getOption("cl.cores", ncores))
  old_best_mse <- Inf # Previous estimate of best forecast MSE
  
  Z_train <- Z[1:(nrow(Z) - window_size),]
  
  
  if (missing(ks)){
    odpc_fit <- crit.odpc(Z=Z_train, k_list=k_list, tol=tol, niter_max=niter_max, ncores=ncores)
    k1 <- odpc_fit[[1]]$k1
    k2 <- odpc_fit[[1]]$k2
    k_tot_max <- k1 + k2
    response <- Z_train[(k_tot_max + 1):(nrow(Z_train)),]
    
    sparse_path <- sparse_odpc_priv(Z=Z_train,
                                    resp=response,
                                    num_lambda_in=nlambda,
                                    alpha_en=alpha_en,
                                    a_ini=odpc_fit[[1]]$a,
                                    D_ini=rbind(odpc_fit[[1]]$alpha, odpc_fit[[1]]$B),
                                    k_tot_max=k1+k2,
                                    k1=k1,
                                    k2=k2,
                                    tol=tol,
                                    eps=eps,
                                    niter_max=niter_max,
                                    pass_grid=FALSE,
                                    lambda_grid_in=c(0))
    
    sparse_path <- lapply(sparse_path, function(fit) {convert_rename_comp(fit, wrap=TRUE, sparse=TRUE)})
    # append original fit to the beginning of the path
    odpc_fit[[1]]$lambda <- 0
    sparse_path <- append(sparse_path, list(odpc_fit), after=0)
    sparse_path <- lapply(sparse_path, function(fit) {construct.odpcs(fit, Z, fn_call=match.call())})
    
    # TODO test is.odpc this
    forecasts <- forecast_sparse_rolled(sparse_path=sparse_path, Z=Z, window_size=window_size, h=h, cl=cl)
    best_fit <- get_best_sparse_fit(forecasts=forecasts, sparse_path=sparse_path, Z=Z, h=h)
  
    opt_comp <- best_fit$opt_comp
    new_best_mse <- best_fit$opt_mse
    new_opt_lambda <- opt_comp[[1]]$lambda
    print('lambda')
    print(new_opt_lambda)
    response_full <- Z[(k_tot_max + 1):(nrow(Z)),]
    if(new_opt_lambda > 0){
      final_fit <- sparse_odpc_priv(Z=Z,
                                    resp=response_full,
                                    num_lambda_in=0,
                                    alpha_en=alpha_en,
                                    a_ini=opt_comp[[1]]$a,
                                    D_ini=rbind(opt_comp[[1]]$alpha, opt_comp[[1]]$B),
                                    k_tot_max=k1+k2,
                                    k1=k1,
                                    k2=k2,
                                    tol=tol,
                                    eps=eps,
                                    niter_max=niter_max,
                                    pass_grid=TRUE,
                                    lambda_grid_in=c(new_opt_lambda))
      final_fit <- convert_rename_comp(final_fit[[1]], wrap=TRUE, sparse=TRUE)
      final_fit <- construct.odpcs(final_fit, data=Z, fn_call=match.call())
    } else {
      final_fit <- crit.odpc(Z=Z, k_list=opt_comp[[1]]$k1, tol=tol, niter_max=niter_max, ncores=ncores)
    }
  } else {
    odpc_fit <- odpc(Z=Z_train, ks=ks, tol=tol, niter_max=niter_max)
    k1 <- odpc_fit[[1]]$k1
    k2 <- odpc_fit[[1]]$k2
    k_tot_max <- k1 + k2
    response <- Z_train[(k_tot_max + 1):(nrow(Z_train)),]
    
    sparse_path <- sparse_odpc_priv(Z=Z_train,
                                    resp=response,
                                    num_lambda_in=nlambda,
                                    alpha_en=alpha_en,
                                    a_ini=odpc_fit[[1]]$a,
                                    D_ini=rbind(odpc_fit[[1]]$alpha, odpc_fit[[1]]$B),
                                    k_tot_max=k1+k2,
                                    k1=k1,
                                    k2=k2,
                                    tol=tol,
                                    eps=eps,
                                    niter_max=niter_max,
                                    pass_grid=FALSE,
                                    lambda_grid_in=c(0))
    
    sparse_path <- lapply(sparse_path, function(fit) {convert_rename_comp(fit, wrap=TRUE, sparse=TRUE)})
    # append original fit to the beginning of the path
    odpc_fit[[1]]$lambda <- 0
    sparse_path <- append(sparse_path, list(odpc_fit), after=0)
    sparse_path <- lapply(sparse_path, function(fit) {construct.odpcs(fit, Z, fn_call=match.call())})
    
    # TODO test is.odpc this
    forecasts <- forecast_sparse_rolled(sparse_path=sparse_path, Z=Z, window_size=window_size, h=h, cl=cl)
    best_fit <- get_best_sparse_fit(forecasts=forecasts, sparse_path=sparse_path, Z=Z, h=h)
    
    opt_comp <- best_fit$opt_comp
    new_best_mse <- best_fit$opt_mse
    new_opt_lambda <- opt_comp[[1]]$lambda
    print('lambda')
    print(new_opt_lambda)
    response_full <- Z[(k_tot_max + 1):(nrow(Z)),]
    if(new_opt_lambda > 0){
      final_fit <- sparse_odpc_priv(Z=Z,
                                    resp=response_full,
                                    num_lambda_in=0,
                                    alpha_en=alpha_en,
                                    a_ini=opt_comp[[1]]$a,
                                    D_ini=rbind(opt_comp[[1]]$alpha, opt_comp[[1]]$B),
                                    k_tot_max=k1+k2,
                                    k1=k1,
                                    k2=k2,
                                    tol=tol,
                                    eps=eps,
                                    niter_max=niter_max,
                                    pass_grid=TRUE,
                                    lambda_grid_in=c(new_opt_lambda))
      final_fit <- convert_rename_comp(final_fit[[1]], wrap=TRUE, sparse=TRUE)
      final_fit <- construct.odpcs(final_fit, data=Z, fn_call=match.call())
    } else {
      final_fit <- odpc(Z=Z, ks=ks, tol=tol, niter_max=niter_max)
    }
  }
  
  on.exit(stopCluster(cl))
  return(final_fit)
}

get_best_sparse_fit <- function(sparse_path, forecasts, Z, h){
  mses <- get_ave_mse_sparse(forecasts=forecasts, Z=Z, h=h)
  print(mses)
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
  fores_comps <- apply(new_comps, 2, function(x, h) { auto <- auto.arima(x, stationary=TRUE, seasonal=FALSE, approximation=TRUE)
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