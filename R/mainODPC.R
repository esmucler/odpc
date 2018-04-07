odpc <- function(Z, ks, method, ini = 'classic', tol = 1e-04, niter_max = 500) {
  # Computes One Sided Dynamic Principal Components.
  #INPUT
  # Z: data matrix, series by columns 
  # ks:  Matrix or vector of integers. If matrix, each row represents the number of lags to use
  # for each component. First column of each row is the number of lags used to define the
  # dynamic principal component (k1), second column is the number of lags of the dynamic principal
  # component used to reconstruct the series (k2). If ks is a vector, its i-th entry 
  # is taken as both k1 and k2 for the i-th component
  # method: Algorithm used. Options are 'ALS' or 'mix'. See details below.
  # ini: initial estimator for the iterations. Either 'classic' or 'gdpc'. Default is 'classic'.
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
    if (!(method %in% c('ALS', 'mix'))) {
      stop("method should be ALS or mix")
    }
  }
  
  num_comp <- nrow(ks)
  
  if (missing(method)){
    if (ncol(Z) > 40){
      method <- 2 # Use mix method for moderately fat data sets
    } else {
      method <- 1
    }
  } else {
    method <- switch(method, 'ALS' = 1, 'mix' = 2)  
  }
  
  output <- vector('list', num_comp)
  k_tot <- apply(ks, 1, sum)
  cum_k_tot_max <- cummax(k_tot) #cumulative maxes of the first reconstructible periods
  k_tot_max_old <- cum_k_tot_max[1] #current max
  res <- Z[((k_tot_max_old)+1):nrow(Z),]  #residuals for the null model
  
  if (ini == 'classic'){
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
  } else if (ini == 'gdpc'){
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
      fit_gdpc <- gdpc(Z = resp, k = ks[iter, 1], tol = tol, niter_max = niter_max)
      f_ini <- fit_gdpc$f[(k_tot_max - ks[iter, 2] + 1):length(fit_gdpc$f)]
      output[[iter]] <- convert_rename_comp(odpc_priv(Z = Z, resp=resp, k_tot_max=k_tot_max,
                                                      k1 = ks[iter, 1], k2 = ks[iter, 2], num_comp=iter, f_ini = f_ini,
                                                      passf_ini = TRUE, tol = tol, niter_max = niter_max, method = method),
                                            wrap=FALSE)
      res <- output[[iter]]$res #current residuals
    }
  }
  criters_conv <- sapply(output, function(x) { return(x$criter) })
  if (any(criters_conv < -1e-5)){
    warning('The last iteration did not decrease the reconstruction MSE')
  }
  fn_call <- match.call()
  output <- construct.odpcs(output, Z, fn_call)
  return(output)
}


#' @title Automatic choosing of tuning parameters for odpc via cross-validation
#' @param Z Data matrix. Each column is a different time series.
#' @param h Forecast horizon.
#' @param k_list List of values of k to choose from.
#' @param max_num_comp Maximum possible number of components to compute
#' @param window_size The size of the rolling window used to estimate the forecasting error.
#' @param ncores_k Number of cores to parallelise over the list of ks.
#' @param ncores_w Number of cores to parallelise over the rolling window (nested in the ks).
#' @param method A string specifying the algorithm used. Options are 'ALS' or 'mix'. See details below.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.
#' @param train_tol Relative precision used in cross-validation. Default is 1e-2.
#' @param train_niter_max Integer. Maximum number of iterations used in cross-validation. Default is 100.


#' @return An object of class odpcs, that is, a list of length equal to the number of computed components, each computed using the optimal value of k. 
#' 
#' @description
#' Computes One-Sided Dynamic Principal Components, choosing the number of components and lags automatically, to minimize an estimate
#' of the forecasting mean squared error.

#' @seealso \code{\link{odpc}}, \code{\link{forecast.odpcs}}
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
    if (!(method %in% c('ALS', 'mix'))) {
      stop("method should be ALS or mix")
    }
  }
  
  if (missing(method)){
    if (ncol(Z) > 40){
      method <- 2 # Use mix method for moderately fat data sets
    } else {
      method <- 1
    }
  } else {
    method <- switch(method, 'ALS' = 1, 'mix' = 2)  
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
      updated_k_params <- update_k_params(ks=ks, k_list=k_list)
      k_maxs <- updated_k_params$k_maxs
      k_trun <- updated_k_params$k_trun
      current_k_max <- updated_k_params$current_k_max
    } else {
      old_best_mse <- new_best_mse
    }
  }
  methods <- c('ALS', 'mix')
  method <- methods[method]
  output <- odpc(Z=Z, ks=as.numeric(ks), method=method, ini = "classic", tol=tol, niter_max=niter_max)
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


#' @title Automatic choosing of tuning parameters for odpc via the minimization of an information criterion
#' @param Z Data matrix. Each column is a different time series.
#' @param k_list List of values of k to choose from.
#' @param max_num_comp Maximum possible number of components to compute
#' @param ncores Number of cores to use in parallel computations.
#' @param method A string specifying the algorithm used. Options are 'ALS' or 'mix'. See details below.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.

#' @return An object of class odpcs, that is, a list of length equal to the number of computed components, each computed using the optimal value of k. 
#' 
#' @description
#' Computes One-Sided Dynamic Principal Components, choosing the number of components and lags automatically, to minimize an 
#' information criterion

#' @seealso \code{\link{odpc}}, \code{\link{cv.odpc}}, \code{\link{forecast.odpcs}}
crit.odpc <- function(Z, k_list = 1:5, max_num_comp = 5, ncores, method, tol = 1e-04, niter_max = 500) {
  
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
    if (!(method %in% c('ALS', 'mix'))) {
      stop("method should be ALS or mix")
    }
  }
  
  if (missing(method)){
    if (ncol(Z) > 40){
      method <- 2 # Use mix method for moderately fat data sets
    } else {
      method <- 1
    }
  } else {
    method <- switch(method, 'ALS' = 1, 'mix' = 2)  
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
  best_fit <- get_best_fit_crit(fits, num_comp=1)
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
    best_fit <- get_best_fit_crit(fits, num_comp=num_comp + 1)
    
    cand_opt_comp <- best_fit$opt_comp
    
    if (best_fit$opt_crit < new_best_crit){
      old_best_crit <- new_best_crit
      new_best_crit <- best_fit$opt_crit
      new_opt_k <- k_list[best_fit$opt_ind]
      res <- opt_comp[[1]]$res
      ks <- c(ks, new_opt_k)
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
