odpc <- function(Z, ks, method, tol = 1e-04, niter_max = 500) {
  # Computes One Sided Dynamic Principal Components.
  #INPUT
  # Z: data matrix, series by columns 
  # ks:  Matrix or vector of integers. If matrix, each row represents the number of lags to use
  # for each component. First column of each row is the number of lags used to define the
  # dynamic principal component (k1), second column is the number of lags of the dynamic principal
  # component used to reconstruct the series (k2). If ks is a vector, its i-th entry 
  # is taken as both k1 and k2 for the i-th component
  # method: Algorithm used. Options are 'ALS', 'mix' or 'gradient'. See details below.
  # niter_max : maximum number of iterations. Default is 500
  # tol : tolerance parameter to declare convergence. Default is 1e-4
  
  #OUTPUT
  # A list of length equal to the number of computed components. The i-th entry of this list is an object of class
  # odpc, that is, a list with entries:
  # a: vector to construct the i-th principal component
  # f: Coordinates of the i-th Principal Component corresponding to the periods k+1,…,T
  # B: matrix of loadings correspoding to f
  # alpha: vector of intercepts corresponding to f
  # mse: mean (in T and m) squared error of the residuals of the fit with the first i components
  # ks: same as input
  # call: the matched call
  # conv: Logical. Did the iterations converge?
  
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
  if (all(!inherits(ks, "matrix"), !inherits(ks, "numeric"))) {
    stop("ks should be a matrix or a vector")
  } else if (!inherits(ks, "matrix")) {
      ks <- cbind(ks, ks)
  }
  if (any(!(ks == floor(ks)), ks < 0)) {
    stop("ks should only have non-negative integers as entries")
  } else if (any(apply(ks, 1, function(x){all(x==0)}))){
    stop("For each component, either k1 or k2 should be positive")
  }
  if (!missing(method)) {
    if (!(method %in% c('ALS', 'mix', 'gradient'))) {
      stop("method should be ALS, mix or gradient")
    }
  }
  
  num_comp <- nrow(ks)
  
  if (missing(method)){
    if (ncol(Z) > 20){
      method <- 3 # Use gradient method for moderately fat data sets
    } else {
      method <- 1
    }
  } else {
    method <- switch(method, 'ALS' = 1, 'mix' = 2, 'gradient' = 3)  
  }
  
  output <- vector('list', num_comp)
  k_tot <- apply(ks, 1, sum)
  cum_k_tot_max <- cummax(k_tot) #cumulative maxes of the first reconstructible periods
  k_tot_max_old <- cum_k_tot_max[1] #current max
  res <- Z[((k_tot_max_old)+1):nrow(Z),]  #residuals for the null model
  
  
  for (iter in seq_len(num_comp)){
    k_tot_max <- cum_k_tot_max[iter]
    # if the max at this iter is larger than the current max
    # we have to lose the difference in periods; else the 
    # reconstructible periods remain the same
    if (k_tot_max > k_tot_max_old){
      resp <- res[((k_tot_max - k_tot_max_old) + 1):nrow(res),]  #current response
      k_tot_max_old <- k_tot_max
    } else {
      resp <- res #current response
    }
    output[[iter]] <- convert_rename_comp(odpc_priv(Z = Z, resp=resp, k_tot_max=k_tot_max,
                                                    k1 = ks[iter, 1], k2 = ks[iter, 2], num_comp=iter, f_ini = 0,
                                                    passf_ini = FALSE, tol = tol, niter_max = niter_max, method = method),
                                          wrap=FALSE)
    res <- output[[iter]]$res #current residuals
  }
  
  criters_conv <- sapply(output, function(x) { return(x$criter) })
  if (any(criters_conv < -1e-5)){
    warning('The last iteration did not decrease the reconstruction MSE')
  }
  fn_call <- match.call()
  output <- construct.odpcs(output, Z, fn_call)
  return(output)
}


#' @title Automatic Choice of Tuning Parameters for One-Sided Dynamic Principal Components via Cross-Validation
#' @param Z Data matrix. Each column is a different time series.
#' @param h Forecast horizon.
#' @param k_list List of values of k to choose from.
#' @param max_num_comp Maximum possible number of components to compute.
#' @param window_size The size of the rolling window used to estimate the forecasting error.
#' @param ncores_k Number of cores to parallelise over \code{k_list}.
#' @param ncores_w Number of cores to parallelise over the rolling window (nested in \code{k_list}).
#' @param method A string specifying the algorithm used. Options are 'ALS', 'mix' or 'gradient'. See details in \code{\link{odpc}}.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.
#' @param train_tol Relative precision used in cross-validation. Default is 1e-2.
#' @param train_niter_max Integer. Maximum number of iterations used in cross-validation. Default is 100.


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
#'\code{components}, \code{fitted}, \code{plot} and \code{print} methods are available for this class.
#' 
#' @description
#' Computes One-Sided Dynamic Principal Components, choosing the number of components and lags automatically, to minimize an estimate
#' of the forecasting mean squared error.
#'
#' @details 
#' We assume that for each component
#' \eqn{k_{1}^{i}=k_{2}^{i}}, that is, the number of lags of \eqn{\mathbf{z}_{t}} used to
#' define the dynamic principal component and the number of lags of
#' \eqn{\widehat{f}^{i}_{t}} used to reconstruct the original series are the same. The number of components and lags
#' is chosen to minimize the cross-validated forecasting error in a
#' stepwise fashion.  
#' Suppose we want to make \eqn{h}-steps ahead forecasts.
#' Let \eqn{w=} \code{window_size}.
#' Then given \eqn{k\in} \code{k_list} we compute the first ODPC
#' defined using \eqn{k} lags, using periods \eqn{1,\dots,T-h-t+1} for \eqn{t=1,\dots,w}, and for each
#' of these fits we compute an h-steps ahead forecast and the corresponding
#' mean squared error \eqn{E_{t,h}}. The cross-validation estimate of the forecasting error
#' is then
#'\deqn{
#'  \widehat{MSE}_{1,k}=\frac{1}{w}\sum\limits_{t=1}^{w}E_{t,h}.
#'}
#'  We choose for the first component the value \eqn{k^{\ast,1}} that minimizes \eqn{\widehat{MSE}_{1,k}}.
#' Then, we fix the first component computed with \eqn{k^{\ast,1}} lags and repeat the
#' procedure with the second component. If the optimal cross-validated
#' forecasting error using the two components, \eqn{\widehat{MSE}_{2,k^{\ast,2}}} is larger than the one using only
#' one component, \eqn{\widehat{MSE}_{1,k^{\ast,1}}}, we stop and output as a final model the ODPC computed using one component
#' defined with \eqn{k^{\ast,1}} lags; otherwise, if \code{max_num_comp} \eqn{\geq 2} we add the second component defined using \eqn{k^{\ast,2}} lags and proceed as before.
#'
#' This method can be computationally costly, especially for large values of  the \code{window_size}. Ideally, the user should set
#' \code{n_cores_k} equal to the length of \code{k_list} and \code{n_cores_w} equal to \code{window_size}; this would entail using
#' \code{n_cores_k} times \code{n_cores_w} cores in total.
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
#' fit <- cv.odpc(x, h=1, k_list = 1:2, max_num_comp = 1, window_size = 2, ncores_k = 2, ncores_w = 2)

cv.odpc <- function(Z, h, k_list = 1:5, max_num_comp = 5, window_size, ncores_k=1, ncores_w=1, method, tol = 1e-04, niter_max = 500, train_tol = 1e-2, train_niter_max = 100) {
  
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
    if (ncol(Z) > 20){
      method <- 3 # Use gradient method for moderately fat data sets
    } else {
      method <- 1
    }
  } else {
    method <- switch(method, 'ALS' = 1, 'mix' = 2, 'gradient' = 3)  
  }
  if (missing(window_size)){
    window_size <- floor(0.2 * nrow(Z))
  }
  
  ks <- c() # List of estimated optimal k for each component
  old_best_mse <- Inf # Previous estimate of best forecast MSE
  
  usable_cores <- min(detectCores(), ncores_k)
  cl <- makeCluster(usable_cores)
  registerDoParallel(cl)
  num_comp <- 1
  
  # one entry per window, roll the training data
  data_field <- build_data_field(Z=Z, window_size = window_size, h = h)
  # one entry per k, obtained by truncating the first 2k peridos of each element in data_field
  response_field <- build_response_field(data_field=data_field, k_trun = 2 * k_list)
  fits <- grid_odpc(data_field = data_field, response_field = response_field,
                    num_comp = num_comp, k_list=k_list, k_maxs = 2 * k_list, window_size=window_size, tol=train_tol,
                    niter_max=train_niter_max, method=method, ncores_w=ncores_w)
  best_fit <- get_best_fit(fits, Z=Z, h=h, window_size = window_size)
  opt_comp <- best_fit$opt_comp
  new_best_mse <- best_fit$opt_mse
  new_opt_k <- k_list[best_fit$opt_ind]
  ks <- c(ks, new_opt_k)
  updated_k_params <- update_k_params(ks=ks, k_list=k_list)
  k_maxs <- updated_k_params$k_maxs
  k_trun <- updated_k_params$k_trun
  current_k_max <- updated_k_params$current_k_max
  
  while (new_best_mse < old_best_mse & num_comp < max_num_comp){
    # response field is now formed by the residuals from the fit using the current optimal componentss
    residual_field <- build_data_field(opt_comp)
    response_field <- build_response_field(data_field=residual_field, k_trun=k_trun)
    # compute another component using the previous fitted ones
    fits <- grid_odpc(data_field = data_field, response_field = response_field,
                      num_comp = num_comp + 1, k_list=k_list, k_maxs=k_maxs, window_size=window_size, tol=train_tol,
                      niter_max=train_niter_max, method=method, ncores_w=ncores_w)
    # append to current components the new fitted ones; extended fits has one entry per k; each of these
    # entries has window_size entries; in each of these we have: the current optimal component computed along the
    # rolling window, with the latest component appended at the end
    extended_fits <- new_window_object(fits, opt_comp)

    # get the optimal k for the new component
    best_fit <- get_best_fit(extended_fits, Z=Z, h=h, window_size = window_size)

    opt_comp <- best_fit$opt_comp

    if (best_fit$opt_mse < new_best_mse){
      old_best_mse <- new_best_mse
      new_best_mse <- best_fit$opt_mse
      new_opt_k <- k_list[best_fit$opt_ind]
      ks <- c(ks, new_opt_k)
      num_comp <- length(ks)
      updated_k_params <- update_k_params(ks=ks, k_list=k_list)
      k_maxs <- updated_k_params$k_maxs
      k_trun <- updated_k_params$k_trun
      current_k_max <- updated_k_params$current_k_max
    } else {
      old_best_mse <- new_best_mse
    }
  }
  methods <- c('ALS', 'mix', 'gradient')
  method <- methods[method]
  output <- odpc(Z=Z, ks=as.numeric(ks), method=method, tol=tol, niter_max=niter_max)
  on.exit(stopCluster(cl))
  return(output)
}

grid_odpc <- function(data_field, response_field, k_list, k_maxs, num_comp, window_size, tol, niter_max, method, ncores_w=1){
    output <- list()
    ind <- NULL
    output <- foreach(ind=1:length(k_list), .packages=c('odpc'))%dopar%{    
                  roll_odpc(data_field=data_field, response_field=response_field[[ind]],
                            k=k_list[ind], k_tot_max=k_maxs[ind], num_comp=num_comp, window_size=window_size, tol=tol,
                            niter_max=niter_max, method=method, ncores=ncores_w)
    }
    output <- convert_rename(output)
    return(output)
}


#' @title Automatic Choice of Tuning Parameters for One-Sided Dynamic Principal Components via the Minimization of an Information Criterion
#' @param Z Data matrix. Each column is a different time series.
#' @param k_list List of values of k to choose from.
#' @param max_num_comp Maximum possible number of components to compute.
#' @param ncores Number of cores to use in parallel computations.
#' @param method A string specifying the algorithm used. Options are 'ALS', 'mix' or 'gradient'. See details in \code{\link{odpc}}.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.

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
#'\code{components}, \code{fitted}, \code{plot} and \code{print} methods are available for this class.

#'
#' @description
#' Computes One-Sided Dynamic Principal Components, choosing the number of components and lags automatically, to minimize an 
#' information criterion.
#' 
#' @details 
#' 
#' We apply the same stepwise approach taken in \code{\link{cv.odpc}}, but now to minimize an
#' information criterion instead of the cross-validated forecasting error. The criterion is
#' inspired by the \eqn{IC_{p3}} criterion proposed in Bai and Ng (2002).
#' Let \eqn{\widehat{\sigma}^{2}_{1,k}} be the reconstruction mean squared error for
#' the first ODPC defined using \eqn{k} lags. Let \eqn{T^{\ast,1,k}=T-2k}.
#' Then we choose the
#' value \eqn{k^{\ast,1}} in \code{k_list} that minimizes
#' \deqn{
#'  {BNG}_{1,k}=\log\left( \widehat{\sigma}^{2}_{1,k} \right)
#'  + ( k+1 ) \frac{\log\left(\min(T^{\ast,1,k},m)\right)}{\min(T^{\ast,1,k},m)}.
#'  }
#' Suppose now that \code{max_num_comp} \eqn{\geq 2} and we
#' have computed \eqn{q-1} dynamic principal components, \eqn{q-1 <} \code{max_num_comp}, each with \eqn{k_{1}^{i}=k_{2}^{i}=k^{\ast, i}} lags, \eqn{i=1,\dots,q-1}.
#' Let \eqn{\widehat{\sigma}^{2}_{q,k}} be the reconstruction mean squared error for
#' the fit obtained using \eqn{q} components, where the first \eqn{q-1} components are defined using 
#' \eqn{k^{\ast, i}}, \eqn{i=1,\dots,q-1} and the last component is defined using \eqn{k} lags.
#' Let \eqn{T^{\ast,q,k}=T-\max\lbrace 2k^{\ast,1},\dots,2k^{\ast,q-1},2k \rbrace}.
#' Let \eqn{k^{\ast,q}} be the value in \code{k_list} that minimizes
#' \deqn{
#'  {BNG}_{q,k}=\log\left( \widehat{\sigma}^{2}_{q,k} \right)
#'  + \left(\sum_{i=1}^{q-1}(k^{\ast,i}+1) + k+1 \right) \frac{\log\left(\min(T^{\ast,q,k},m)\right)}{\min(T^{\ast,q,k},m)}  .
#'  }
#' If \eqn{{BNG}_{q,k^{\ast,q}}} is larger than \eqn{{BNG}_{q-1,k^{\ast,q-1}}}
#' we stop and the final model is the ODPC with \eqn{q-1} components. Else we add the \eqn{q}-th component defined using \eqn{k^{\ast,q}}
#' and continue as before.
#' 
#' @references
#' Peña D., Smucler E. and Yohai V.J. (2017). “Forecasting Multiple Time Series with One-Sided Dynamic Principal Components.” Available at https://arxiv.org/abs/1708.04705.
#' 
#' Bai J. and Ng S. (2002). “Determining the Number of Factors in Approximate Factor Models.” Econometrica, 70(1), 191–221.

#' @seealso \code{\link{odpc}}, \code{\link{cv.odpc}}, \code{\link{forecast.odpcs}}
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
#' fit <- crit.odpc(x, k_list = 1:2, max_num_comp = 1)

crit.odpc <- function(Z, k_list = 1:5, max_num_comp = 5, ncores = 1, method, tol = 1e-04, niter_max = 500) {
  
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
    if (ncol(Z) > 20){
      method <- 3 # Use gradient method for moderately fat data sets
    } else {
      method <- 1
    }
  } else {
    method <- switch(method, 'ALS' = 1, 'mix' = 2, 'gradient' = 3)  
  }
  
  ks <- c() # List of estimated optimal k for each component
  old_best_crit <- Inf # Previous estimate of best forecast MSE
  
  usable_cores <- min(detectCores(), ncores)
  cl <- makeCluster(usable_cores)
  registerDoParallel(cl)
  num_comp <- 1
  k_trun <- 2 * k_list
  
  response_field <- lapply(k_trun, function(k){ Z[(k+1):nrow(Z),]  })
  fits <- grid_crit_odpc(Z = Z, response_field = response_field, k_maxs = 2 * k_list, num_comp = num_comp,
                         k_list=k_list, tol=tol, niter_max=niter_max, method=method)
  best_fit <- get_best_fit_crit(fits, k_acum = 0)
  opt_comp <- best_fit$opt_comp
  new_best_crit <- best_fit$opt_crit
  new_opt_k <- k_list[best_fit$opt_ind]
  res <- opt_comp[[1]]$res
  ks <- c(ks, new_opt_k)
  updated_k_params <- update_k_params(ks=ks, k_list=k_list)
  k_maxs <- updated_k_params$k_maxs
  k_trun <- updated_k_params$k_trun
  current_k_max <- updated_k_params$current_k_max

  
  while (new_best_crit < old_best_crit & num_comp < max_num_comp){
    
    # compute another component using the previous fitted ones
    response_field <- lapply(k_trun, function(k){ res[(k+1):nrow(res),]  })
    fits <- grid_crit_odpc(Z = Z, response_field = response_field, k_maxs = k_maxs, num_comp = num_comp + 1,
                           k_list=k_list, tol=tol, niter_max=niter_max, method=method)
    
    
    # get the optimal k for the new component
    best_fit <- get_best_fit_crit(fits, k_acum = sum(ks + 1))
    
    cand_opt_comp <- best_fit$opt_comp
    
    if (best_fit$opt_crit < new_best_crit){
      old_best_crit <- new_best_crit
      new_best_crit <- best_fit$opt_crit
      new_opt_k <- k_list[best_fit$opt_ind]
      opt_comp <- c(opt_comp, cand_opt_comp)
      res <- cand_opt_comp[[1]]$res
      ks <- c(ks, new_opt_k)
      num_comp <- length(ks)
      updated_k_params <- update_k_params(ks=ks, k_list=k_list)
      k_maxs <- updated_k_params$k_maxs
      k_trun <- updated_k_params$k_trun
      current_k_max <- updated_k_params$current_k_max
    } else {
      old_best_crit <- new_best_crit
    }
  }
  
  methods <- c('ALS', 'mix')
  method <- methods[method]
  fn_call <- match.call()
  output <- construct.odpcs(opt_comp, Z, fn_call)
  on.exit(stopCluster(cl))
  return(output)
}

grid_crit_odpc <- function(Z, response_field, k_maxs, k_list, num_comp, tol, niter_max, method){
  output <- list()
  ind <- NULL
  output <- foreach(ind=1:length(k_list), .packages=c('odpc'))%dopar%{    
      convert_rename_comp(odpc_priv(Z = Z, resp = response_field[[ind]], k_tot_max = k_maxs[ind],
                                    k1 = k_list[ind], k2 = k_list[ind], num_comp=num_comp,
                                    f_ini = c(0), passf_ini = FALSE, tol = tol, niter_max = niter_max,
                                    method = method),
                          wrap=TRUE)
  }
  return(output)
}
