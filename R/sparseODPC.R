#' @title Automatic Choice of Regularization Parameters for Sparse One-Sided Dynamic Principal Components using a BIC type criterion
#' @param Z Data matrix. Each column is a different time series.
#' @param k_list List of values of k to choose from.
#' @param max_num_comp Maximum possible number of components to compute.
#' @param nlambda Length of penalty sequence.
#' @param tol Relative precision. Default is 1e-4.
#' @param niter_max Integer. Maximum number of iterations. Default is 500.
#' @param eps Between 0 and 1, used to build penalty sequence
#' @param ncores Number of cores to use in parallel computations



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
#'\item{lambda}{Regularization parameter used for this component.}
#'\code{components}, \code{fitted}, \code{plot} and \code{print} methods are available for this class.
#' 
#' 
#' @description
#' Computes Sparse One-Sided Dynamic Principal Components, choosing the number of components and regularization parameters automatically, using a BIC type criterion.
#'
#' @details 
#' First \code{\link{crit.odpc}} is called to choose the number of lags and of components to use. Each sparse component is then computed using a regularized version of the
#' odpc objective function (see \code{\link{odpc}}), where the L1 norm of the \eqn{\mathbf{a}} vector is penalized. The penalization parameter \eqn{\lambda} is chosen from a grid of candidates
#' of size \code{nlambda}, seeking to minimize the following BIC type criterion
#' \deqn{
#' \log(MSE(\mathbf{a}_{\lambda},\mathbf{\alpha}_{\lambda}, \mathbf{B}_{\lambda} )) + \frac{\log(T^{\ast} m)}{T^{\ast}m} \Vert \mathbf{a}_{\lambda}\Vert_{0},
#' }
#' where \eqn{\mathbf{a}_{\lambda},\mathbf{B}_{\lambda} } are the estimates associated with a given \eqn{\lambda}, \eqn{m} is the number of series and
#' \eqn{T^{\ast}} is the number of periods being reconstructed.
#' 
#' @references
#' Peña D., Smucler E. and Yohai V.J. (2019). “Forecasting Multiple Time Series with One-Sided Dynamic Principal Components.” Journal of the American Statistical Association.
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
#' fit <- crit.sparse_odpc(x, k_list = 1, ncores = 1)
crit.sparse_odpc <- function(Z,  k_list = 1:3, max_num_comp=1, nlambda=20, tol = 1e-04, niter_max = 500, eps=1e-3, ncores=1) {
  
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
  
  odpc_fit <- crit.odpc(Z=Z, max_num_comp=max_num_comp, k_list=k_list, tol=tol, niter_max=niter_max, ncores=ncores)
  num_comp_total <- length(odpc_fit)
  k1 <- odpc_fit[[1]]$k1
  k2 <- odpc_fit[[1]]$k2
  k_tot_max <- k1 + k2
  response <- Z[(k_tot_max + 1):(nrow(Z)),]
  num_comp <- 1
  fit_res <- get_partial_comp(Z=Z, response=response, 
                              nlambda=nlambda, k1=k1, k2=k2, k_tot_max=k_tot_max, tol=tol,
                              eps=eps, niter_max=niter_max, 
                              odpc_fit=list(odpc_fit[[num_comp]]))
  final_fit <- fit_res$final_fit
  while(num_comp < num_comp_total){
    num_comp <- num_comp + 1
    k1 <- odpc_fit[[num_comp]]$k1
    k2 <- odpc_fit[[num_comp]]$k2
    k_tot_max <- max(k1 + k2, k_tot_max)
    response <- Z[(k_tot_max + 1):(nrow(Z)),] - fitted(final_fit, num_comp=num_comp)
    fit_res <- get_partial_comp(Z=Z, response=response, 
                                nlambda=nlambda, k1=k1, k2=k2, k_tot_max=k_tot_max, tol=tol,
                                eps=eps, niter_max=niter_max, 
                                odpc_fit=list(odpc_fit[[num_comp]]))
    final_fit <- append(final_fit, fit_res$final_fit)
  }
    
  return(final_fit)
}


get_partial_comp <- function(Z, response, nlambda, k_tot_max,
                             k1, k2, tol, eps, niter_max, odpc_fit){
  
  sparse_path <- sparse_odpc_priv(Z=Z,
                                  resp=response,
                                  num_lambda_in=nlambda,
                                  a_ini=odpc_fit[[1]]$a,
                                  D_ini=rbind(odpc_fit[[1]]$alpha, odpc_fit[[1]]$B),
                                  k_tot_max=k_tot_max,
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
  
  mses <- sapply(sparse_path, function(x) x[[1]]$mse)
  numb_vars <- sapply(sparse_path, function(x) sum(x[[1]]$a != 0))
  Tast <- nrow(response)
  m <- ncol(response)
  bics <- log(mses) + log(Tast * m) / (Tast * m) * numb_vars
  
  opt_ind <- which.min(bics)
  opt_comp <- sparse_path[[opt_ind]]
  
  
  
  return(list(final_fit=opt_comp))
}