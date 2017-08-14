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
  res <- Z
  
  if (ini == 'classic'){
    for (iter in seq_len(num_comp)){
      output[[iter]] <- odpc_priv(Z = res,  k1 = ks[iter, 1], k2 = ks[iter, 2], f_ini = 0,
                                  passf_ini = FALSE, tol = tol, niter_max = niter_max,
                                  method = method)
      res <- output[[iter]]$res
    }
  } else if (ini == 'gdpc'){
    for (iter in seq_len(num_comp)){
      fit_gdpc <- gdpc(Z = res, k = ks[iter, 1], tol = tol, niter_max = niter_max)
      f_ini <- fit_gdpc$f[(ks[iter, 1] + 1):length(fit_gdpc$f)]
      # f_ini <- fit_gdpc$f[1:(length(fit_gdpc$f) - ks[iter, 1])]
      output[[iter]] <- odpc_priv(Z = res, k1 = ks[iter, 1], k2 = ks[iter, 2], f_ini = f_ini,
                                   passf_ini = TRUE, tol = tol, niter_max = niter_max,
                                   method = method)
      res <- output[[iter]]$res
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